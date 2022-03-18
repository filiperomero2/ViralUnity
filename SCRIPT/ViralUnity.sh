#!/bin/bash

## ViralUnity.sh v.1.0.2-beta ##

# Written by Moreira and D'arc 2022.

##############################################################################################################################
# Set functions

help(){

echo "
ViralUnity

Description: ViralUnity is a pipeline for the inference of viral consensus genome sequences from Illumina paired end reads.

It takes as arguments:

	--LIBDIR - Absolute path for libraries' root directory;
	--OUTDIR - Absolute path for output directory;
	--REF - Path for a reference genome in fasta format;
	--ADAPTERS - Absolute path for trimmomatic adapters fasta file;
	--MINCOV - Minimum sequencing coverage to call a base on the consensus sequence (Optional; default = 100);
	--MINLEN - Minimum read length (Optional; default = 50);
	--HEADCROP - Number of bases to trim from the start of the read, useful for primer sequences removal (Optional; default = 30);
	--THREADS - Number of threads available for processing (Optional; default = 1).

Minimal usage:$ ./ViralUnity.sh --LIBDIR ~/LIBRARIES/RUN_1/ --OUTDIR ~/ANALYSIS/RUN1/ --REF ~/REFERENCE_GENOMES/reference.fasta --ADAPTERS ~/trimmomatic/adapter.fa

Alternatively: $ ./ViralUnity.sh --LIBDIR ~/LIBRARIES/RUN_1/ --OUTDIR ~/ANALYSIS/RUN1/ --REF ~/REFERENCE_GENOMES/reference.fasta --ADAPTERS ~/trimmomatic/adapter.fa --MINCOV 200 --MINLEN 40 --HEADCROP 20 --THREADS 6
"
}

version(){
	echo "ViralUnity v1.0.2-beta"
}

clean(){
	rm $NAME.temp.consensus.fa
	mv trim.* $OUTDIR/FILTERED_DATA/
	mv *sorted.bam $OUTDIR/MAPPING/
	rm *bam
	mv *bai $OUTDIR/MAPPING/
	mv $NAME.calls* $OUTDIR/MAPPING/
	mv *.txt  $OUTDIR/MAPPING/
	mv masked.$NAME.consensus.fa $OUTDIR/CONSENSUS/
}

writestats(){
	# Number of masked bases
	NBASES=$(grep -v '>' masked.$NAME.consensus.fa | tr -cd 'N' | wc -c)
	# Genome size
	LEN=$(wc -l $NAME.table_cov_basewise.txt | awk '{print $1}')
	# Genome coverage
	GENCOV=$(calc 1-$NBASES/$LEN)
	# Number of raw reads
	RAW=$(grep -cE "^\+$" *fastq | head -n 1 | sed -E 's/.+\://g')
	RAW=$(calc $RAW*2)
	# Number of paired filtered reads
	PAIRED=$(grep -cE "^\+$" trim.p*fastq | sed -E 's/.+\://g' | head -n 1)
	PAIRED=$(calc $PAIRED*2)
	# Number of unpaired filtered reads
	UNPAIRED=$(grep -cE "^\+$" $U | sed -E 's/.+\://g')
	# Number of reads mapped
	MAPPED=$(samtools view -c -F 260 $NAME.sorted.bam)
	# Relative frequency of mapped reads
	USAGE=$(calc $MAPPED/$RAW)
	# Average depth
	SUM=$(cat  $NAME.table_cov_basewise.txt | awk '{sum+=$3; print sum}' | tail -n 1)
	COV=$(calc $SUM/$LEN)
	# Percentage of genome coverage at 10x
	TEMP=$(awk  '$3 > 10' $NAME.table_cov_basewise.txt | wc -l)
	COV10=$(calc $TEMP/$LEN)
	# Percentage of genome coverage at 100x
	TEMP=$(awk  '$3 > 100' $NAME.table_cov_basewise.txt | wc -l)
	COV100=$(calc $TEMP/$LEN)
	# Percentage of genome coverage at 1000x
	TEMP=$(awk  '$3 > 1000' $NAME.table_cov_basewise.txt | wc -l)
	COV1000=$(calc $TEMP/$LEN)
	# Print statistics to .csv file
	echo "$NAME,$RAW,$PAIRED,$UNPAIRED,$MAPPED,$USAGE,$COV,$COV10,$COV100,$COV1000,$GENCOV" >> ../$FILE.stats_report.csv

}

calc() { awk "BEGIN{print $*}"; }


##############################################################################################################################
# Validate arguments

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	help
	exit
fi

if [ -z "$1" ] || [[ $1 == -v ]] || [[ $1 == --version ]]; then
	version
	exit
fi

echo ""
echo "###############################"
echo "Step 1: Arguments Verification"
echo "###############################"
echo ""


while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        echo $1 $2
   fi
   shift
done


if [[ -d $LIBDIR ]] 
then 
	echo "Root directory for samples indentified ->  $LIBDIR"
else 
	echo "Samples directory not found. Please specify."
	exit
fi


if [[ -d $OUTDIR ]] 
then 
	echo "Output directory indentified ->  $OUTDIR"
	echo "Please specify another to avoid overwriting previous analysis."
	exit
else 
	echo "Samples directory created -> $OUTDIR"
	mkdir $OUTDIR
fi


if [[ -f $REF ]] 
then 
	echo "Reference genome file indentified ->  $REF"
else 
	echo "Reference genome file not found. Please specify."
	exit
fi

if [[ -f $ADAPTERS ]] 
then 
	echo "Trimmomatic adapters fasta file identified ->  $ADAPTERS"
else 
	echo "Trimmomatic adapters fasta file not found. Please specify."
	exit
fi

if [ -z "$MINCOV" ];
then
      echo "Minimum coverage is not defined, using the default (100x)"
	  MINCOV=100;
else
      echo "Minimum coverage threshold identified: $MINCOV"
fi

if [ -z "$MINLEN" ]
then
      echo "Minimum read length not specified, using the default (50bp)"
	  MINLEN=50
else
      echo "Minimum read length: $MINLEN"
fi

if [ -z "$HEADCROP" ]
then
      echo "Number of bases to trim from the beggining of read not specified, using the default (30bp)"
	  HEADCROP=30
else
      echo "Number of bases to trim from the beggining of read: $HEADCROP"
fi


if [ -z "$THREADS" ]
then
      echo "Number of threads not specified, starting with the default (1)"
	  THREADS=1
else
      echo "Number of threads available for processing: $THREADS"
fi


##############################################################################################################################
# Start the pipeline

echo ""
echo "###############################"
echo "Step 2: Start Assembly"
echo "###############################"
echo ""

echo ""
date


# Index reference genome 
FILE=$(echo $REF | sed 's/.*\///g')
echo "Reference file name -> $FILE"
REFNAME=$(cat $REF | head -n 1 | sed 's/>//' | sed 's/\s.*//')
REFPATH=$(echo $REF | sed "s/$FILE//g")

cd $REFPATH
bowtie2-build -q $REF reference
REF2=$(pwd)/reference


# Go to libraries' root directory and set $OUTDIR/*/
cd $LIBDIR

date > $FILE.timereport.txt

mkdir -p $OUTDIR/QC_RAW/ $OUTDIR/QC_FILTERED/ $OUTDIR/FILTERED_DATA/ $OUTDIR/MAPPING/ $OUTDIR/CONSENSUS/

# Print the header of sequencing_stats.csv
echo "sample,number_of_raw_reads,number_of_paired_filtered_reads,number_of_unpaired_filtered_reads,number_of_mapped_reads,efficiency,average_depth,coverage_10x,coverage_100x,coverage_1000x,genome_coverage" > $FILE.stats_report.csv

# For each sample
for SAMPLE in */
do
	# Skip $OUTDIR/
	if [ $SAMPLE == "$OUTDIR/" ]
	then
		continue
	fi

	# Enter sample diretory
	cd $SAMPLE

	# Print sample directory name
	echo ""
	echo "Working on -> $SAMPLE"
    date
	echo "Working on -> $SAMPLE" >> ../$FILE.timereport.txt
	date >> ../$FILE.timereport.txt
	echo ""
			


	# Decompress data
	gunzip *gz
	
	# Check raw data files
	N=$(ls *fastq | wc -l)
	if [ $N == 2 ]
	then
		R1="$(ls *fastq | head -1)"
		R2="$(ls *fastq | tail -1)"
	else
		echo "Please provide only two fastq files (R1 and R2) for each sample."
		echo "Only files with the .fastq extension are used."
		echo "Different number of files found in directory $SAMPLE."
		exit
	fi
	
	# Get sample name. It should never include a _ character (used as separator)
	NAME=$(echo $R1 | sed -E 's/_.+//g')

	# QC report for raw data 
	fastqc -q -t $THREADS *fastq
	mv *zip $OUTDIR/QC_RAW/
	mv *html $OUTDIR/QC_RAW/
		
	# Filter data with fastp
	echo "Performing strict QC..."
	if [ $HEADCROP -eq 0 ]
	then
		trimmomatic PE -threads $THREADS -phred33 $R1 $R2 trim.p.$R1 trim.u.$R1 trim.p.$R2 trim.u.$R2 ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:$MINLEN
	else
		trimmomatic PE -threads $THREADS -phred33 $R1 $R2 trim.p.$R1 trim.u.$R1 trim.p.$R2 trim.u.$R2 ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 HEADCROP:$HEADCROP MINLEN:$MINLEN
	fi

	# QC report for filtered data
	fastqc -q -t $THREADS trim*fastq
	mv *html $OUTDIR/QC_FILTERED/
	mv *zip $OUTDIR/QC_FILTERED/
	
	# Concatenate unpaired reads
	cat trim.u.*fastq > trim.uni.$R1

	# Get the names of qc reads
	U=trim.uni.$R1
	R1=trim.p.$R1
	R2=trim.p.$R2

	# Map reads against reference genome and handle the bam file
	echo "Mapping..."
	bowtie2 --very-sensitive --no-unal -p $THREADS -x $REF2 -1 $R1 -2 $R2 -U $U | samtools view -bS - > $NAME.bam   
    samtools faidx $REF 
	samtools sort $NAME.bam -o $NAME.sorted.bam
	samtools index $NAME.sorted.bam 
	
	# Call variants 
	echo "Calling variants..."
	bcftools mpileup --threads $THREADS --max-depth 20000 -E -Ou -f $REF $NAME.sorted.bam  | bcftools call --ploidy 1 --threads $THREADS -mv -Oz -o $NAME.calls.vcf.gz
	bcftools norm --threads $THREADS -f $REF $NAME.calls.vcf.gz -Oz -o $NAME.calls.norm.vcf.gz
	bcftools filter --threads $THREADS --IndelGap 5 $NAME.calls.norm.vcf.gz -Oz -o $NAME.calls.norm.flt-indels.vcf.gz
    bcftools index --threads $THREADS $NAME.calls.norm.flt-indels.vcf.gz
	
	# Get consensus sequences
	echo "Inferring consensus sequences..."
	bcftools consensus -f $REF $NAME.calls.norm.flt-indels.vcf.gz > $NAME.temp.consensus.fa
	bedtools genomecov -bga -ibam  $NAME.sorted.bam > $NAME.table_cov.txt
	bedtools genomecov -d -ibam  $NAME.sorted.bam > $NAME.table_cov_basewise.txt
	awk  '$4 < '$MINCOV'' $NAME.table_cov.txt > $NAME.table_coverage_min-cov-$MINCOV.txt
	bedtools maskfasta -fi  $NAME.temp.consensus.fa -fo masked.$NAME.consensus.fa -bed $NAME.table_coverage_min-cov-$MINCOV.txt
	sed -i "s/$REFNAME/$NAME/"  masked.$NAME.consensus.fa

	# Compute statistics and write to the stats report
	writestats

    # Compress fastq and move files to their respective directories
	gzip *fastq
	clean

	echo ""
	date
	date  >> ../$FILE.timereport.txt
	
	# Go to the previous directory
	cd ../

done

##############################################################################################################################
# Main loop finished, setting final results. 

mv $FILE.stats_report.csv $OUTDIR/
mv $FILE.timereport.txt $OUTDIR/

cd $OUTDIR/

echo ""
echo "###############################"
echo "Step 3: Final Compilation"
echo "###############################"
echo ""

# Concatenate consensus sequences
cat CONSENSUS/*fa > CONSENSUS/CONSENSUS.fasta

# Multiqc plot for raw data
multiqc QC_RAW/
mv multiqc_report.html multiqc_report_RAW.html
mv multiqc_* QC_RAW/

# Multiqc plot for filtered data
multiqc QC_FILTERED/
mv multiqc_report.html multiqc_report_FILTERED.html
mv multiqc_* QC_FILTERED/



echo ""
echo "#######################"
echo "####### The end #######"
echo "#######################"
echo ""

date
echo "Finished" >> $FILE.timereport.txt
date  >> $FILE.timereport.txt
echo ""

cd

exit
