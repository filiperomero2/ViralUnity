#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
from typing import Dict

import pandas as pd
import yaml


def count_fastq_reads(path: str) -> int:
    """
    Count reads in FASTQ/FASTQ.GZ by counting lines and dividing by 4.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"FASTQ not found: {path}")

    opener = gzip.open if path.endswith(".gz") else open
    n_lines = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for _ in fh:
            n_lines += 1

    if n_lines == 0:
        return 0

    if n_lines % 4 != 0:
        sys.stderr.write(
            f"WARNING: FASTQ line count not divisible by 4: {path} (lines={n_lines}).\n"
        )

    return n_lines // 4


def load_samples_from_config(config_path: str) -> Dict[str, str]:
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    if "samples" not in cfg:
        raise ValueError(f"Config file '{config_path}' must contain a 'samples' key.")

    if not isinstance(cfg["samples"], dict):
        raise ValueError("'samples' in config must be a dict: sample -> FASTQ path.")

    return cfg["samples"]


def add_rpm(
    df: pd.DataFrame,
    sample_to_fastq: Dict[str, str],
    reads_col: str,
    rpm_col: str,
) -> pd.DataFrame:

    if "sample" not in df.columns:
        raise ValueError("Input summary TSV must contain a 'sample' column.")
    if reads_col not in df.columns:
        raise ValueError(
            f"Input summary TSV missing reads column '{reads_col}'. "
            f"Available columns: {list(df.columns)}"
        )

    totals = {}
    for s in sorted(df["sample"].unique()):
        if s not in sample_to_fastq:
            raise KeyError(f"Sample '{s}' not found in config['samples'].")
        totals[s] = count_fastq_reads(sample_to_fastq[s])

    out = df.copy()
    out["total_reads"] = out["sample"].map(totals)

    denom = out["total_reads"].astype(float)
    numer = out[reads_col].astype(float)

    out[rpm_col] = 0.0
    mask = denom > 0
    out.loc[mask, rpm_col] = (numer.loc[mask] / denom.loc[mask]) * 1e6

    return out


def run_cli():
    ap = argparse.ArgumentParser(
        description="Add RPM to ViralUnity taxa summary TSV using raw reads as denominator."
    )
    ap.add_argument("--summary", required=True, help="Input taxa summary TSV.")
    ap.add_argument("--out", required=True, help="Output TSV with RPM column added.")
    ap.add_argument("--config", required=True, help="ViralUnity config YAML.")
    ap.add_argument(
        "--reads-col",
        default="mapped_reads",
        help="Numerator column (e.g. mapped_reads for contigs, count for reads).",
    )
    ap.add_argument(
        "--rpm-col",
        default="rpm",
        help="RPM column name to create (default: rpm).",
    )
    args = ap.parse_args()

    sample_to_fastq = load_samples_from_config(args.config)
    df = pd.read_csv(args.summary, sep="\t")

    out = add_rpm(
        df,
        sample_to_fastq,
        reads_col=args.reads_col,
        rpm_col=args.rpm_col,
    )
    out.to_csv(args.out, sep="\t", index=False)


def run_snakemake():
    in_path = str(snakemake.input[0])
    out_path = str(snakemake.output[0])

    sample_to_fastq = getattr(snakemake.params, "samples", None)
    if sample_to_fastq is None:
        raise ValueError(
            "Snakemake params must include 'samples=config[\"samples\"]'."
        )

    reads_col = getattr(snakemake.params, "reads_col", "mapped_reads")
    rpm_col = getattr(snakemake.params, "rpm_col", "rpm")

    df = pd.read_csv(in_path, sep="\t")
    out = add_rpm(
        df,
        sample_to_fastq,
        reads_col=reads_col,
        rpm_col=rpm_col,
    )
    out.to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
