#!/usr/bin/env python3
"""
summarize_krona_taxa.py

Standardized taxonomic summary from Krona input files (seq_id<TAB>taxid).

Inputs:
  - One Krona input file (2-column TSV: sequence_id, taxid)
  - NCBI taxdump directory (nodes.dmp, names.dmp)

Outputs:
  - taxa_summary.tsv with schema:
      sample, tool, mode, rank, taxid, name, count, percent, source

Notes:
  - tool: kraken2 | diamond
  - mode: reads | contigs
  - rank: family | genus | species
  - percent is computed as count / total_assigned * 100

Supports:
  - CLI execution (argparse)
  - Snakemake `script:` execution via the `snakemake` object
"""

import argparse
import os
import csv
from collections import defaultdict


############################
# Taxdump utilities
############################

def load_taxdump(taxdump_dir: str):
    nodes_path = os.path.join(taxdump_dir, "nodes.dmp")
    names_path = os.path.join(taxdump_dir, "names.dmp")

    if not os.path.isfile(nodes_path):
        raise FileNotFoundError(f"nodes.dmp not found at: {nodes_path}")
    if not os.path.isfile(names_path):
        raise FileNotFoundError(f"names.dmp not found at: {names_path}")

    parent = {}
    rank = {}
    with open(nodes_path, newline="") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            taxid = parts[0]
            parent[taxid] = parts[1]
            rank[taxid] = parts[2]

    name = {}
    with open(names_path, newline="") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            taxid = parts[0]
            nm = parts[1]
            cls = parts[3]
            if cls == "scientific name":
                name[taxid] = nm

    return parent, rank, name


def lineage_at_ranks(taxid: str, parent, rank_map, name_map, wanted):
    """
    Walk up taxonomy until root, collecting first hit for wanted ranks.
    Returns dict rank -> (taxid, name)
    """
    out = {r: (None, None) for r in wanted}
    cur = taxid
    visited = set()

    while cur and cur not in visited:
        visited.add(cur)
        r = rank_map.get(cur)
        if r in wanted and out[r][0] is None:
            out[r] = (cur, name_map.get(cur, ""))
        cur = parent.get(cur)

    return out


############################
# Core summarization logic
############################

def summarize_krona_file(
    krona_input: str,
    taxdump_dir: str,
    tool: str,
    mode: str,
    sample: str,
    output_tsv: str,
):
    parent, rank_map, name_map = load_taxdump(taxdump_dir)

    wanted_ranks = ["family", "genus", "species"]
    counts = defaultdict(int)
    total = 0

    # Fail-safe: if input is missing/empty, write header-only output and return
    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)

    if (not os.path.exists(krona_input)) or (os.path.getsize(krona_input) == 0):
        with open(output_tsv, "w", newline="") as out:
            out.write("sample\ttool\tmode\trank\ttaxid\tname\tcount\tpercent\tsource\n")
        return

    with open(krona_input, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            taxid = str(row[1]).strip()
            if not taxid or taxid == "0":
                continue

            total += 1
            lineage = lineage_at_ranks(taxid, parent, rank_map, name_map, wanted_ranks)
            for r, (tid, _) in lineage.items():
                if tid:
                    counts[(r, tid)] += 1

    # IMPORTANT: write a fresh per-sample output (not append)
    with open(output_tsv, "w", newline="") as out:
        out.write("sample\ttool\tmode\trank\ttaxid\tname\tcount\tpercent\tsource\n")
        for (r, tid), cnt in sorted(counts.items()):
            nm = name_map.get(tid, "")
            pct = (cnt / total * 100.0) if total > 0 else 0.0
            out.write(
                f"{sample}\t{tool}\t{mode}\t{r}\t{tid}\t{nm}\t{cnt}\t{pct:.4f}\t{krona_input}\n"
            )


############################
# Entrypoints
############################

def main_cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxdump", required=True, help="Path to NCBI taxdump directory")
    ap.add_argument("--tool", required=True, choices=["kraken2", "diamond"])
    ap.add_argument("--mode", required=True, choices=["reads", "contigs"])
    ap.add_argument("--sample", required=True, help="Sample ID")
    ap.add_argument("--input", required=True, help="Krona input TSV (seq_id<TAB>taxid)")
    ap.add_argument("--output", required=True, help="Output taxa_summary.tsv")
    args = ap.parse_args()

    summarize_krona_file(
        krona_input=args.input,
        taxdump_dir=args.taxdump,
        tool=args.tool,
        mode=args.mode,
        sample=args.sample,
        output_tsv=args.output,
    )


def main_snakemake():
    # Snakemake injects a global `snakemake` object
    krona_input = str(snakemake.input[0])
    output_tsv = str(snakemake.output[0])

    taxdump_dir = str(snakemake.params["taxdump"])
    tool = str(snakemake.params["tool"])
    mode = str(snakemake.params["mode"])

    # sample can be passed explicitly, or derived from wildcards
    sample = str(snakemake.params.get("sample", snakemake.wildcards.sample))

    summarize_krona_file(
        krona_input=krona_input,
        taxdump_dir=taxdump_dir,
        tool=tool,
        mode=mode,
        sample=sample,
        output_tsv=output_tsv,
    )


if __name__ == "__main__":
    # If executed via Snakemake `script:`, a global `snakemake` exists
    if "snakemake" in globals():
        main_snakemake()
    else:
        main_cli()
