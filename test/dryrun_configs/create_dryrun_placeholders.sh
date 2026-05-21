#!/bin/bash

# Create dummy directories
mkdir -p data/reads
mkdir -p data/references
mkdir -p databases/kraken2
mkdir -p databases/krona/taxonomy
mkdir -p databases/taxdump
mkdir -p databases/diamond

# Create dummy FASTQ files
touch data/reads/SAMPLE1_R1.fastq.gz
touch data/reads/SAMPLE1_R2.fastq.gz
touch data/reads/SAMPLE2_R1.fastq.gz
touch data/reads/SAMPLE2_R2.fastq.gz
touch data/reads/SAMPLE_NP.fastq.gz
touch data/reads/SAMPLE_NP_SEG.fastq.gz
touch data/reads/META1_R1.fastq.gz
touch data/reads/META1_R2.fastq.gz
touch data/reads/META_NP.fastq.gz

# Create dummy reference files
touch data/references/reference.fasta
touch data/references/segment1.fasta
touch data/references/segment2.fasta
touch data/references/segment3.fasta

# Create dummy database files
touch databases/kraken2/hash.k2d
touch databases/krona/taxonomy/taxonomy.tab
touch databases/taxdump/nodes.dmp
touch databases/taxdump/names.dmp
touch databases/diamond/nr.dmnd

echo "Placeholder files created successfully in data/, databases/, etc."
echo "You can now run 'snakemake -n' for any of the dry-run configs."

mkdir -p databases/virus_genomes
touch databases/virus_genomes/viral.genomes.fasta
touch databases/virus_genomes/genome2taxid.tsv
