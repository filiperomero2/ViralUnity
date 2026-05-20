#!/usr/bin/env python

import gzip
import subprocess

import pandas as pd


def get_number_of_reads(fastq: str) -> int:
    """Count records in a FASTQ file by counting the "+" separator lines.

    Works for both plain and gzip-compressed inputs. Avoids subprocess +
    shell=True so the path is not interpreted by a shell.
    """
    opener = gzip.open if fastq.endswith(".gz") else open
    count = 0
    with opener(fastq, "rt") as f:
        for line in f:
            if line.rstrip("\n") == "+":
                count += 1
    return count


def get_number_of_mapped_reads(bam: str) -> int:
    """Return the number of mapped reads in a BAM via `samtools view -c -F 260`."""
    result = subprocess.run(
        ["samtools", "view", "-c", "-F", "260", bam],
        check=True,
        capture_output=True,
        text=True,
    )
    return int(result.stdout.strip())


def get_coverage_info(table_cov, minimum_depth):
    df = pd.read_csv(table_cov, header=None, sep=r"\s+")
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


def generate_output(
    fastq, trim_fastq, bam, table_cov, sample_name, minimum_depth, output, segment=None
):
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
    data = [sample_name]
    if segment is not None:
        data.append(segment)

    data.extend(
        [
            number_of_reads,
            number_of_trim_reads,
            number_of_mapped_reads,
            average_depth,
            percentage_of_sites_above_10x,
            percentage_of_sites_above_100x,
            percentage_of_sites_above_1000x,
            percentage_of_sites_above_specified_threshold,
        ]
    )
    return pd.DataFrame(data).T


def main(input, output, minimum_depth, segment=None):
    # keep fist read and ignore paired in illumina scenario
    fastq = input[0]
    trim_fastq = input[2]
    bam = input[3]
    table_cov = input[4]

    sample_name = bam.replace(".sorted.bam", "")
    sample_name = sample_name.split("/")[-1]

    df_out = generate_output(
        fastq, trim_fastq, bam, table_cov, sample_name, minimum_depth, output, segment
    )

    df_out.to_csv(output, header=False, index=False)


if __name__ == "__main__":
    input = snakemake.input
    output = snakemake.output[0]
    minimum_depth = snakemake.params[0]
    segment = getattr(snakemake.wildcards, "segment", None)

    main(input, output, minimum_depth, segment)

    exit()
