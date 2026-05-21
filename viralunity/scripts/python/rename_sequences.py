#!/usr/bin/env python
"""Rename a per-sample consensus FASTA's header to the sample name.

This script is executed by Snakemake via a ``script:`` directive (see
``rules/consensus_*.smk::rename_sequences``). Snakemake injects a
``snakemake`` global at runtime (``snakemake.input``, ``snakemake.output``)
which ruff/mypy cannot see; F821/ignore_errors are configured for this
directory in ``pyproject.toml``.

The script reads a single-sample consensus FASTA (where the original header
may be anything), strips it, and writes a single record whose header is
``>{sample_name}`` derived from the input filename (``.consensus.fasta``
suffix stripped).
"""

import re


def rename_sequences(input_path: str, output: str) -> None:
    """Rewrite the FASTA at ``input_path`` so its single record's header is the sample name.

    The sample name is the basename of ``input_path`` minus the
    ``.consensus.fasta`` suffix. The sequence body is concatenated from
    every non-header line in the source file (newlines preserved as in the
    input).
    """
    sample_name = re.sub(".+/", "", input_path)
    sample_name = re.sub(".consensus.fasta", "", sample_name)
    with open(input_path) as f:
        lines = f.readlines()

    sequence = ""

    for line in lines:
        if not line.startswith(">"):
            sequence = sequence + line

    with open(output, "w") as f:
        renamed_sequence = ">" + sample_name + "\n" + sequence
        f.write(renamed_sequence)


if __name__ == "__main__":
    print("Running rename_sequences.py")
    rename_sequences(snakemake.input[0], snakemake.output[0])
    print("Finished rename_sequences.py")
    exit()
