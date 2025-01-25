#!/usr/bin/env python

import re
import sys
import pandas as pd
import subprocess

input = snakemake.input
fastq = input[0]
trim_fastq = input[2]
bam = input[3]
table_cov = input[4]

output = snakemake.output[0]

minimum_depth = snakemake.params[0]

sample_name = bam.replace(".sorted.bam","")
sample_name = sample_name.split("/")[-1]


def get_number_of_reads(fastq):
    if fastq.endswith('.gz'):
        command = "gunzip -c " + fastq + " | grep -cE \"^\+$\""
    else:
        command = "grep -cE \"^\+$\" " + fastq 
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    number_of_reads = int(out)
    return(number_of_reads)

def get_number_of_mapped_reads(bam):
    command = "samtools view -c -F 260 " + bam
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    number_of_mapped_reads = int(out)
    return(number_of_mapped_reads)

def get_coverage_info(table_cov):
    df =  pd.read_csv(table_cov,header=None,delim_whitespace=True) # check this    
    total_sequenced_bases = df[2].sum()
    reference_genome_length = len(df)
    average_depth = total_sequenced_bases/reference_genome_length
    percentage_of_sites_above_10x = len(df[df[2]>=10])/reference_genome_length
    percentage_of_sites_above_100x = len(df[df[2]>=100])/reference_genome_length
    percentage_of_sites_above_1000x = len(df[df[2]>=1000])/reference_genome_length
    percentage_of_sites_above_specified_threshold = len(df[df[2]>=minimum_depth])/reference_genome_length
    return(average_depth,percentage_of_sites_above_10x,percentage_of_sites_above_100x,percentage_of_sites_above_1000x,percentage_of_sites_above_specified_threshold)


def write_output():
    number_of_reads = get_number_of_reads(fastq)
    number_of_trim_reads = get_number_of_reads(trim_fastq)
    number_of_mapped_reads = get_number_of_mapped_reads(bam)
    average_depth,percentage_of_sites_above_10x,percentage_of_sites_above_100x,percentage_of_sites_above_1000x,percentage_of_sites_above_specified_threshold = get_coverage_info(table_cov)
    data = [sample_name,
            number_of_reads,
            number_of_trim_reads,
            number_of_mapped_reads,
            average_depth,
            percentage_of_sites_above_10x,
            percentage_of_sites_above_100x,
            percentage_of_sites_above_1000x,
            percentage_of_sites_above_specified_threshold]
    df = pd.DataFrame(data).T
    df.to_csv(output,header=False,index=False)

write_output()

exit()