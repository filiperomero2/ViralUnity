#!/bin/bash

## GetConSeqVir.sh v.1 ##

# It is a simple pipeline for inference of viral consensus genome sequences from Illumina paired end data.
# The pipeline performs QC, mapping, variant calling and consensus sequence inference using diverse tools.
# It takes 2 positional arguments:
#	1 - Absolute path for libraries' root directory;
#	2 - Path for a reference genome in fasta format;
#	4 - Minimum sequencing coverage to call a base on the consensus sequence (default = 100x)
#	5 - Number of threads available for processing (default = 1)
# Example usage:$ ./ViralUnit.sh ~/LIBRARIES/RUN_1/ ~/REFERENCE_GENOMES/reference.fasta 200 6

##############################################################################################################################

# Set variables...
LIBDIR=$1
REF=$2
MINCOV=$3
THREADS=$4

##############################################################################################################################
# Validate arguments and define a function

if [[ -d $LIBDIR ]] 
then 
	echo "Root directory for samples indentified ->  $LIBDIR"
else 
	echo "Samples directory not found. Please specify."
	exit
fi

if [[ -f $REF ]] 
then 
	echo "Reference genome file indentified ->  $REF"
else 
	echo "Reference genome file not found. Please specify."
	exit
fi

if [ -z "$MINCOV" ];
then
      echo "Minimum coverage is not defined, using the default (100x)"
	  MINCOV=100;
else
      echo "Minimum coverage threshold identified: $MINCOV"
fi

if [ -z "$THREADS" ]
then
      echo "Number of threads not specified, starting with the default (1)"
	  THREADS=1
else
      echo "Number of threads available for processing: $THREADS"
fi

calc() { awk "BEGIN{print $*}"; }

##############################################################################################################################
# Start the pipeline

# Index reference genome 
FILE=$(echo $REF | sed 's/.*\///g')
REFNAME=$(cat $REF | head -n 1 | sed 's/>//')
REFPATH=$(echo $REF | sed "s/$FILE//g")
cd $REFPATH
bowtie2-build $REF reference
REF2=$(pwd)/reference


# Go to libraries' root directory and set $FILE.RESULTS/*/
cd $LIBDIR
mkdir -p $FILE.RESULTS/RAW_QC/ $FILE.RESULTS/FILTERED_QC/ $FILE.RESULTS/FILTERED_DATA/ $FILE.RESULTS/MAPPING/ $FILE.RESULTS/CONSENSUS/

# Print the header of sequencing_stats.csv
echo "sample,number_of_raw_reads,number_of_paired_filtered_reads,number_of_unpaired_filtered_reads,number_of_mapped_reads,efficiency,average_depth,coverage_10x,coverage_100x,coverage_1000x,genome_coverage" > stats_report.csv

# For each sample
for SAMPLE in */
do
	# Skip $FILE.RESULTS/
	if [ $SAMPLE == "$FILE.RESULTS/" ]
	then
		continue
	fi

	# Print sample directory name
	echo "$SAMPLE"
		
	# Enter sample diretory
	cd $SAMPLE

	# Decompress data
	gunzip *gz
	
	# Check raw data files
	N=$(ls  | wc -l)
	if [ $N == 2 ]
	then
		R1="$(ls | head -1)"
		R2="$(ls | tail -1)"
	else
		echo "Please provide only two fastq files (R1 and R2) for each sample."
		echo "Different number of files found in directory $SAMPLE"
		exit
	fi
	
	NAME=$(echo $R1 | sed -E 's/_.+//g')

	# QC report for raw data 
	fastqc -t $THREADS *fastq
	mv *zip ../$FILE.RESULTS/RAW_QC/
	mv *html ../$FILE.RESULTS/RAW_QC/
		
	# Filter data with fastp
	echo "Performing strict QC..."
	fastp --detect_adapter_for_pe --thread $THREADS -i $R1 -I $R2 -o trim.p.$R1 -O trim.p.$R2 --unpaired1 trim.u.$R1 --unpaired2 trim.u.$R2
	mv fastp.html $NAME.fastp.html
	mv fastp.json $NAME.fastp.json

	# QC report for filtered data
	fastqc -t $THREADS trim*fastq
	mv *json ../$FILE.RESULTS/FILTERED_QC/
	mv *html ../$FILE.RESULTS/FILTERED_QC/
	mv *zip ../$FILE.RESULTS/FILTERED_QC/
	
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
	bcftools mpileup --max-depth 20000 -E -Ou -f $REF $NAME.sorted.bam  | bcftools call --ploidy 1 -mv -Oz -o $NAME.calls.vcf.gz
	bcftools norm -f $REF $NAME.calls.vcf.gz -Oz -o $NAME.calls.norm.vcf.gz
	bcftools filter --IndelGap 5 $NAME.calls.norm.vcf.gz -Oz -o $NAME.calls.norm.flt-indels.vcf.gz
    bcftools index $NAME.calls.norm.flt-indels.vcf.gz
	
	# Get consensus sequences
	echo "Inferring consensus sequences..."
	cat $REF | bcftools consensus $NAME.calls.norm.flt-indels.vcf.gz > $NAME.temp.consensus.fa
	bedtools genomecov -bga -ibam  $NAME.sorted.bam > $NAME.table_cov.txt
	bedtools genomecov -d -ibam  $NAME.sorted.bam > $NAME.table_cov_basewise.txt
	awk  '$4 < '$MINCOV'' $NAME.table_cov.txt > $NAME.table_coverage_min-cov-$MINCOV.txt
	bedtools maskfasta -fi  $NAME.temp.consensus.fa -fo masked.$NAME.consensus.fa -bed $NAME.table_coverage_min-cov-$MINCOV.txt
	sed -i "s/$REFNAME/$NAME/"  masked.$NAME.consensus.fa

	# Compute statistics and write to the stats report
	NBASES=$(grep -v '>' masked.$NAME.consensus.fa | tr -cd 'N' | wc -c)
	LEN=$(wc -l $NAME.table_cov_basewise.txt | awk '{print $1}')
	GENCOV=$(calc 1-$NBASES/$LEN)
	RAW=$(grep -cE "^\+$" *fastq | head -n 1 | sed -E 's/.+\://g')
	RAW=$(calc $RAW*2)
	PAIRED=$(grep -cE "^\+$" trim.p*fastq | sed -E 's/.+\://g' | head -n 1)
	PAIRED=$(calc $PAIRED*2)
	UNPAIRED=$(grep -cE "^\+$" $U | sed -E 's/.+\://g')
	MAPPED=$(samtools view -c -F 260 $NAME.sorted.bam)
	USAGE=$(calc $MAPPED/$RAW)
	SUM=$(cat  $NAME.table_cov_basewise.txt | awk '{sum+=$3; print sum}' | tail -n 1)
	COV=$(calc $SUM/$LEN)
	TEMP=$(awk  '$3 > 10' $NAME.table_cov_basewise.txt | wc -l)
	COV10=$(calc $TEMP/$LEN)
	TEMP=$(awk  '$3 > 100' $NAME.table_cov_basewise.txt | wc -l)
	COV100=$(calc $TEMP/$LEN)
	TEMP=$(awk  '$3 > 1000' $NAME.table_cov_basewise.txt | wc -l)
	COV1000=$(calc $TEMP/$LEN)
	echo "$NAME,$RAW,$PAIRED,$UNPAIRED,$MAPPED,$USAGE,$COV,$COV10,$COV100,$COV1000,$GENCOV" >> ../stats_report.csv

    	# Move files to their respective directories
	gzip *fastq
	rm $NAME.temp.consensus.fa
	mv trim.* ../$FILE.RESULTS/FILTERED_DATA/
	mv *bam ../$FILE.RESULTS/MAPPING/
	mv *bai ../$FILE.RESULTS/MAPPING/
	mv $NAME.calls* ../$FILE.RESULTS/MAPPING/
	mv *.txt  ../$FILE.RESULTS/MAPPING/
	mv masked.$NAME.consensus.fa ../$FILE.RESULTS/CONSENSUS/

	echo ""
	
	# Go to the previous directory
	cd ../

done

##############################################################################################################################
# Main loop finished, setting final results. 

mv stats_report.csv $FILE.RESULTS/

cd $FILE.RESULTS/

# Concatenate consensus sequences
cat CONSENSUS/*fa > CONSENSUS/CONSENSUS.fasta

# Multiqc plot for raw data
multiqc RAW_QC/
mv multiqc_report.html multiqc_report_RAW.html
mv multiqc_* RAW_QC/

# Multiqc plot for filtered data
multiqc FILTERED_QC/
mv multiqc_report.html multiqc_report_FILTERED.html
mv multiqc_* FILTERED_QC/

cd ..

echo "#######################"
echo "####### The end #######"
echo "#######################"

exit
