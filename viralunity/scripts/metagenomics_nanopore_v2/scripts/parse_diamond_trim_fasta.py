#!/usr/bin/env python3
"""Trim contig FASTA sequences to viral‑matching regions reported by DIAMOND.

The script can run **stand‑alone** (using regular CLI flags) **or** be executed via
Snakemake’s `script:` directive.  It auto‑detects the environment so you don’t have
to maintain two versions.

Input (either mode):
  1. A contig FASTA file – the same file queried with DIAMOND (nucleotide seqs).
  2. A DIAMOND tabular output (.tsv) that includes *qstart* and *qend* (cols 7‑8).

It gathers all (qstart, qend) ranges per contig, merges overlaps/adjacencies,
extracts the nucleotide subsequences, and writes them to one FASTA where each entry
is named `<contig>|region<N>_<start>_<end>`.

CLI usage example ────────────────────────────────────────────────────────────────
    python parse_diamond_trim_fasta.py \
        --fasta final.contigs.fa \
        --diamond sample.diamond.tsv \
        --output trimmed_contigs.fa \
        --min-len 150

Snakemake rule snippet ───────────────────────────────────────────────────────────
    rule trim_diamond_contig_regions:
        input:
            fasta   = "{sample}.fa",
            tsv     = "{sample}.diamond.tsv"
        output:
            trimmed = "{sample}.diamond.trimmed.fa"
        params:
            min_len = 150
        log: "logs/{sample}.trim.log"
        script: "workflow/scripts/parse_diamond_trim_fasta.py"

The `log:` entry will automatically contain everything printed by this script.
"""


import argparse
from pathlib import Path
from typing import Dict, List, Tuple

################################################################################
# Argument parsing (CLI mode)                                                  #
################################################################################

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Extract viral‑matching regions from contigs based on DIAMOND output.")
    p.add_argument("--fasta", required=True, help="Input contig FASTA (nucleotide)")
    p.add_argument("--diamond", required=True, help="DIAMOND tabular output (.tsv)")
    p.add_argument("--output", required=True, help="Output FASTA with trimmed regions")
    p.add_argument(
        "--min-len",
        type=int,
        default=1,
        help="Discard trimmed regions shorter than this length (default: 1)")
    return p.parse_args()

################################################################################
# FASTA helpers (no external dependencies required)                            #
################################################################################

def read_fasta(path: Path) -> Dict[str, str]:
    """Return a dict {seq_id: sequence}."""
    seqs: Dict[str, List[str]] = {}
    sid: str | None = None
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                sid = line[1:].split()[0]
                seqs[sid] = []
            elif sid:
                seqs[sid].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def write_fasta(records: List[Tuple[str, str]], path: Path, width: int = 60) -> None:
    """Write `(id, seq)` tuples to *path* wrapping to *width* chars/line."""
    with path.open("w") as fh:
        for sid, seq in records:
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")

################################################################################
# DIAMOND helpers                                                              #
################################################################################

def parse_diamond_tsv(path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """Return {qseqid: [(qstart, qend), ...]} from DIAMOND output."""
    hits: Dict[str, List[Tuple[int, int]]] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            qseqid = parts[0]
            try:
                qstart = int(parts[6])
                qend   = int(parts[7])
            except ValueError:
                continue  # malformed coordinates
            if qstart > qend:  # can be reversed on minus strand
                qstart, qend = qend, qstart
            hits.setdefault(qseqid, []).append((qstart, qend))
    return hits


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping/adjacent (start, end) intervals (1‑based, inclusive)."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged: List[List[int]] = [list(intervals[0])]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1][1] = max(last_end, end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]

################################################################################
# Main                                                                         #
################################################################################

def main() -> None:
    """Entry point – works in both CLI and Snakemake contexts."""
    # ── Detect execution context ────────────────────────────────────────────
    if "snakemake" in globals():  # running via Snakemake `script:`
        fasta_path   = Path(snakemake.input.fasta)
        diamond_path = Path(snakemake.input.tsv)
        output_path  = Path(snakemake.output.trimmed)
        min_len      = int(getattr(snakemake.params, "min_len", 1))
    else:                           # running from the command line
        args         = parse_args()
        fasta_path   = Path(args.fasta)
        diamond_path = Path(args.diamond)
        output_path  = Path(args.output)
        min_len      = args.min_len

    # ── Load data ───────────────────────────────────────────────────────────
    contigs = read_fasta(fasta_path)
    hits    = parse_diamond_tsv(diamond_path)

    # ── Extract & write ─────────────────────────────────────────────────────
    records: List[Tuple[str, str]] = []
    for qid, seq in contigs.items():
        if qid not in hits:
            continue
        for idx, (start, end) in enumerate(merge_intervals(hits[qid]), 1):
            if end - start + 1 < min_len:
                continue
            subseq = seq[start - 1 : end]  # 0‑based slice
            records.append((f"{qid}|region{idx}_{start}_{end}", subseq))

    write_fasta(records, output_path)

    print(f"[✓] Wrote {len(records)} trimmed sequences → {output_path}")


if __name__ == "__main__":
    main()
