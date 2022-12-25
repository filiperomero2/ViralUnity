#!/bin/bash

## ViralUnity.sh v.1.0.5 ##

# Written by Moreira et al 2022.

##############################################################################################################################
# Set functions

help(){

echo "
ViralUnity

Description: ViralUnity is a pipeline for the inference of viral consensus genome sequences from Illumina paired end reads.

It takes as arguments:

	--input - Absolute path for libraries' root directory;
	--output - Absolute path for output directory;
	--reference - Path for a reference genome in fasta format;
	--adapters - Absolute path for trimmomatic adapters fasta file;
	--minimum_depth - Minimum sequencing depth to call a base on the consensus sequence (Optional; default = 100);
	--minimum_length - Minimum read length (Optional; default = 50);
	--trim - Number of bases to trim from the start of the read, useful for primer sequences removal (Optional; default = 30);
	--threads - Number of threads available for processing (Optional; default = 1).

Minimal usage: $ ./viralunity.sh --input ~/LIBRARIES/RUN_1/ --output ~/ANALYSIS/RUN1/ --reference ~/REFERENCE_GENOMES/reference.fasta --adapters ~/trimmomatic/adapter.fa

Alternatively: $ ./viralunity.sh --input ~/LIBRARIES/RUN_1/ --output ~/ANALYSIS/RUN1/ --reference ~/REFERENCE_GENOMES/reference.fasta --adapters ~/trimmomatic/adapter.fa --minimum_depth 200 --minimum_length 40 --trim 20 --threads 6
"
}

version(){
	echo "ViralUnity v1.0.5"
}

clean(){
	rm $name.temp.consensus.fa
	mv trim.* $output/filtered_data/
	mv *sorted.bam $output/mapping_and_variant_call_files/
	mv *bai $output/mapping_and_variant_call_files/
	mv $name.calls* $output/mapping_and_variant_call_files/
	mv *.txt  $output/mapping_and_variant_call_files/
	mv masked.$name.consensus.fa $output/consensus/
}

writestats(){
	# Number of masked bases
	nbases=$(grep -v '>' masked.$name.consensus.fa | tr -cd 'N' | wc -c)
	# Genome size
	length=$(wc -l $name.table_cov_basewise.txt | awk '{print $1}')
	# Genome coverage
	coverage=$(calc 1-$nbases/$length)
	# Number of raw reads
	raw=$(grep -cE "^\+$" *fastq | head -n 1 | sed -E 's/.+\://g')
	raw=$(calc $raw*2)
	# Number of paired filtered reads
	paired=$(grep -cE "^\+$" trim.p*fastq | sed -E 's/.+\://g' | head -n 1)
	paired=$(calc $paired*2)
	# Number of unpaired filtered reads
	unpaired=$(grep -cE "^\+$" $U | sed -E 's/.+\://g')
	# Number of reads mapped
	mapped=$(samtools view -c -F 260 $name.sorted.bam)
	# Relative frequency of mapped reads
	usage=$(calc $mapped/$raw)
	# Average depth
	sum=$(cat  $name.table_cov_basewise.txt | awk '{sum+=$3; print sum}' | tail -n 1)
	average_depth=$(calc $sum/$length)
	# Percentage of genome coverage at 10x
	temp=$(awk  '$3 > 10' $name.table_cov_basewise.txt | wc -l)
	cov10=$(calc $temp/$length)
	# Percentage of genome coverage at 100x
	temp=$(awk  '$3 > 100' $name.table_cov_basewise.txt | wc -l)
	cov100=$(calc $temp/$length)
	# Percentage of genome coverage at 1000x
	temp=$(awk  '$3 > 1000' $name.table_cov_basewise.txt | wc -l)
	cov1000=$(calc $temp/$length)
	# Print statistics to .csv file
	echo "$name,$raw,$paired,$Uunpaired,$mapped,$usage,$average_depth,$cov10,$cov100,$cov1000,$coverage" >> ../stats_report.csv

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


if [[ -d $input ]] 
then 
	echo "Root directory for samples indentified ->  $input"
else 
	echo "Samples directory not found. Please specify."
	exit
fi


if [[ -d $output ]] 
then 
	echo "Output directory indentified ->  $output"
	echo "Please specify another to avoid overwriting previous analysis."
	exit
else 
	echo "Samples directory created -> $output"
	mkdir $output
fi


if [[ -f $reference ]] 
then 
	echo "Reference genome file indentified ->  $reference"
else 
	echo "Reference genome file not found. Please specify."
	exit
fi

if [[ -f $adapters ]] 
then 
	echo "Trimmomatic adapters fasta file identified ->  $adapters"
else 
	echo "Trimmomatic adapters fasta file not found. Please specify."
	exit
fi

if [ -z "$minimum_depth" ];
then
      echo "Minimum coverage is not defined, using the default (100x)"
	  minimum_depth=100;
else
      echo "Minimum coverage threshold identified: $minimum_depth"
fi

if [ -z "$minimum_length" ]
then
      echo "Minimum read length not specified, using the default (50bp)"
	  minimum_length=50
else
      echo "Minimum read length: $minimum_length"
fi

if [ -z "$trim" ]
then
      echo "Number of bases to trim from the beggining of read not specified, using the default (30bp)"
	  trim=30
else
      echo "Number of bases to trim from the beggining of read: $trim"
fi


if [ -z "$threads" ]
then
      echo "Number of threads not specified, starting with the default (1)"
	  threads=1
else
      echo "Number of threads available for processing: $threads"
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


# Get reference genome info
reference_file_name=$(echo $reference | sed 's/.*\///g')
echo "Reference file name -> $reference_file_name"
reference_name=$(cat $reference | head -n 1 | sed 's/>//' | sed 's/\s.*//')


# Go to libraries' root directory and set $output/*/
cd $input

date > timereport.txt

mkdir -p $output/qc_report/raw/ $output/qc_report/filtered/ $output/filtered_data/ $output/mapping_and_variant_call_files/ $output/consensus/

# Print the header of sequencing_stats.csv
echo "sample,number_of_raw_reads,number_of_paired_filtered_reads,number_of_unpaired_filtered_reads,number_of_mapped_reads,efficiency,average_depth,coverage_10x,coverage_100x,coverage_1000x,genome_coverage" > stats_report.csv

# For each sample
for sample in */
do
	# Skip $output/
	if [ $sample == "$output/" ]
	then
		continue
	fi

	# Enter sample diretory
	cd $sample

	# Print sample directory name
	echo ""
	echo "Working on -> $sample"
    date
	echo "Working on -> $sample" >> ../timereport.txt
	date >> ../timereport.txt
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
		echo "Different number of files found in directory $sample."
		exit
	fi
	
	# Get sample name. It should never include a _ character (used as separator)
	name=$(echo $R1 | sed -E 's/_.+//g')

	# QC report for raw data 
	fastqc -q -t $threads *fastq
	mv *zip $output/qc_report/raw/
	mv *html $output/qc_report/raw/
		
	# Filter data with fastp
	echo "Performing strict QC..."
	if [ $trim -eq 0 ]
	then
		trimmomatic PE -threads $threads -phred33 $R1 $R2 trim.p.$R1 trim.u.$R1 trim.p.$R2 trim.u.$R2 ILLUMINACLIP:$adapters:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:$minimum_length
	else
		trimmomatic PE -threads $threads -phred33 $R1 $R2 trim.p.$R1 trim.u.$R1 trim.p.$R2 trim.u.$R2 ILLUMINACLIP:$adapters:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 HEADCROP:$trim MINLEN:$minimum_length
	fi

	# QC report for filtered data
	fastqc -q -t $threads trim*fastq
	mv *html $output/qc_report/filtered/
	mv *zip $output/qc_report/filtered/
	
	# Concatenate unpaired reads
	cat trim.u.*fastq > trim.uni.$R1

	# Get the names of qc reads
	U=trim.uni.$R1
	R1=trim.p.$R1
	R2=trim.p.$R2

	# Map reads against reference genome and handle the bam file
	echo "Mapping..."
	minimap2 -a -t $threads -x sr  $reference $R1 $R2 | samtools view -bS -F 4 - | samtools sort -o $name.sorted.bam -
	samtools index $name.sorted.bam 
	
	# Call variants 
	echo "Calling variants..."
	bcftools mpileup --threads $threads --max-depth 20000 -E -Ou -f $reference $name.sorted.bam  | bcftools call --ploidy 1 --threads $threads -mv -Oz -o $name.calls.vcf.gz
	bcftools norm --threads $threads -f $reference $name.calls.vcf.gz -Oz -o $name.calls.norm.vcf.gz
	bcftools filter --threads $threads --IndelGap 5 $name.calls.norm.vcf.gz -Oz -o $name.calls.norm.flt-indels.vcf.gz
    bcftools index --threads $threads $name.calls.norm.flt-indels.vcf.gz
	
	# Get consensus sequences
	echo "Inferring consensus sequences..."
	bcftools consensus -f $reference $name.calls.norm.flt-indels.vcf.gz > $name.temp.consensus.fa
	bedtools genomecov -bga -ibam  $name.sorted.bam > $name.table_cov.txt
	bedtools genomecov -d -ibam  $name.sorted.bam > $name.table_cov_basewise.txt
	awk  '$4 < '$minimum_depth'' $name.table_cov.txt > $name.table_coverage_min-cov-$minimum_depth.txt
	bedtools maskfasta -fi  $name.temp.consensus.fa -fo masked.$name.consensus.fa -bed $name.table_coverage_min-cov-$minimum_depth.txt
	sed -i "s/$reference_name/$name/"  masked.$name.consensus.fa

	# Compute statistics and write to the stats report
	writestats

    # Compress fastq and move files to their respective directories
	gzip *fastq
	clean

	echo ""
	date
	date  >> ../timereport.txt
	
	# Go to the previous directory
	cd ../

done

##############################################################################################################################
# Main loop finished, setting final results. 

mv stats_report.csv $output/
mv timereport.txt $output/

cd $output/

echo ""
echo "###############################"
echo "Step 3: Final Compilation"
echo "###############################"
echo ""

# Concatenate consensus sequences
cat consensus/*fa > consensus/consensus.fasta

# Multiqc plot for raw data
multiqc qc_report/raw/
mv multiqc_report.html multiqc_report_RAW.html
mv multiqc_* qc_report/raw/

# Multiqc plot for filtered data
multiqc qc_report/filtered/
mv multiqc_report.html multiqc_report_FILTERED.html
mv multiqc_* qc_report/filtered/



echo ""
echo "#######################"
echo "####### The end #######"
echo "#######################"
echo ""

date
echo "Finished" >> timereport.txt
date  >> timereport.txt
echo ""

cd

exit
