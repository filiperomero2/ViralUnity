#!/usr/bin/env python
"""Per-sample assembly-statistics summary for the viralunity consensus pipeline.

This script is executed by Snakemake via a ``script:`` directive in
``rules/consensus_*_common.smk``. Snakemake injects a ``snakemake`` global
at runtime (``snakemake.input``, ``snakemake.output``, ``snakemake.params``,
``snakemake.wildcards``) — ruff/mypy cannot see it and are silenced for this
file via the ``viralunity/scripts/python/*.py`` overrides in
``pyproject.toml``.

For each sample the script writes a one-row CSV containing read counts,
mean depth, and the fraction of reference positions exceeding 10x, 100x,
1000x, and the user-supplied ``minimum_depth`` threshold. In segmented mode
the segment name is inserted as the second column.
"""

import gzip
import subprocess
from typing import Optional, Sequence, Tuple

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


def get_coverage_info(
    table_cov: str, minimum_depth: int
) -> Tuple[float, float, float, float, float]:
    """Read a basewise-coverage TSV and return mean depth + fractions of
    positions above 10x / 100x / 1000x / ``minimum_depth``.
    """
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
    fastq: str,
    trim_fastq: str,
    bam: str,
    table_cov: str,
    sample_name: str,
    minimum_depth: int,
    output: str,
    segment: Optional[str] = None,
) -> pd.DataFrame:
    """Assemble the one-row DataFrame written to ``output``."""
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


def main(
    input_paths: Sequence[str],
    output: str,
    minimum_depth: int,
    segment: Optional[str] = None,
) -> None:
    """Snakemake entry point.

    ``input_paths`` is the list passed via ``snakemake.input``; positions 0,
    2, 3, 4 are used (raw fastq, trimmed fastq, sorted BAM, basewise-coverage
    TSV). Position 1 is the Illumina paired R2 fastq, which is intentionally
    ignored — the pipeline reports stats per-sample, not per-strand.
    """
    fastq = input_paths[0]
    trim_fastq = input_paths[2]
    bam = input_paths[3]
    table_cov = input_paths[4]

    sample_name = bam.replace(".sorted.bam", "")
    sample_name = sample_name.split("/")[-1]

    df_out = generate_output(
        fastq, trim_fastq, bam, table_cov, sample_name, minimum_depth, output, segment
    )

    df_out.to_csv(output, header=False, index=False)


if __name__ == "__main__":
    main(
        snakemake.input,
        snakemake.output[0],
        snakemake.params[0],
        getattr(snakemake.wildcards, "segment", None),
    )

    exit()
