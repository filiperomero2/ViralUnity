#!/usr/bin/env python

import pandas as pd
import subprocess


def get_number_of_reads(fastq):
    if fastq.endswith(".gz"):
        command = "gunzip -c " + fastq + ' | grep -cE "^\+$"'
    else:
        command = 'grep -cE "^\+$" ' + fastq
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    number_of_reads = int(out)
    return number_of_reads


def get_number_of_mapped_reads(bam):
    command = "samtools view -c -F 260 " + bam
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    number_of_mapped_reads = int(out)
    return number_of_mapped_reads


def get_coverage_info(table_cov, minimum_depth):
    df = pd.read_csv(table_cov, header=None, delim_whitespace=True)  # check this
    total_sequenced_bases = df[2].sum()
    reference_genome_length = len(df)
    average_depth = total_sequenced_bases / reference_genome_length
    percentage_of_sites_above_10x = len(df[df[2] >= 10]) / reference_genome_length
    percentage_of_sites_above_100x = len(df[df[2] >= 100]) / reference_genome_length
    percentage_of_sites_above_1000x = len(df[df[2] >= 1000]) / reference_genome_length
    percentage_of_sites_above_specified_threshold = (
        len(df[df[2] >= minimum_depth]) / reference_genome_length
    )
    return (
        average_depth,
        percentage_of_sites_above_10x,
        percentage_of_sites_above_100x,
        percentage_of_sites_above_1000x,
        percentage_of_sites_above_specified_threshold,
    )


def generate_output(fastq, trim_fastq, bam, table_cov, sample_name, minimum_depth, output):
    number_of_reads = get_number_of_reads(fastq)
    number_of_trim_reads = get_number_of_reads(trim_fastq)
    number_of_mapped_reads = get_number_of_mapped_reads(bam)
    (
        average_depth,
        percentage_of_sites_above_10x,
        percentage_of_sites_above_100x,
        percentage_of_sites_above_1000x,
        percentage_of_sites_above_specified_threshold,
    ) = get_coverage_info(table_cov, minimum_depth)
    data = [
        sample_name,
        number_of_reads,
        number_of_trim_reads,
        number_of_mapped_reads,
        average_depth,
        percentage_of_sites_above_10x,
        percentage_of_sites_above_100x,
        percentage_of_sites_above_1000x,
        percentage_of_sites_above_specified_threshold,
    ]
    return pd.DataFrame(data).T


def main(input, output, minimum_depth):
    fastq = input[0]
    trim_fastq = input[2]
    bam = input[3]
    table_cov = input[4]

    sample_name = bam.replace(".sorted.bam", "")
    sample_name = sample_name.split("/")[-1]

    df_out = generate_output(fastq, trim_fastq, bam, table_cov, sample_name, minimum_depth, output)
    
    df_out.to_csv(output, header=False, index=False)


if __name__ == "__main__":  # TODO: Check if this brake snakemake
    input = snakemake.input
    output = snakemake.output[0]
    minimum_depth = snakemake.params[0]

    main(input, output, minimum_depth)

    exit()
