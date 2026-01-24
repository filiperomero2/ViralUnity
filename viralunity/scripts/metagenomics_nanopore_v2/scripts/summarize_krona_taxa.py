#!/usr/bin/env python3

import os
import sys
from collections import defaultdict

RANKS_OF_INTEREST = ("family", "genus", "species")


def load_taxdump(nodes_dmp, names_dmp):
    parent = {}
    rank = {}
    name = {}

    with open(nodes_dmp) as f:
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split("|")]

            taxid = parts[0]
            parent_taxid = parts[1]
            rank_name = parts[2]

            parent[taxid] = parent_taxid
            rank[taxid] = rank_name

    with open(names_dmp) as f:
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split("|")]

            taxid = parts[0]
            taxname = parts[1]
            name_class = parts[3]

            if name_class == "scientific name":
                name[taxid] = taxname

    return parent, rank, name



def load_diamond_reads(diamond_tax_file):
    reads = {}
    with open(diamond_tax_file) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            contig_id = parts[0]
            mapped_reads = int(parts[-4])
            reads[contig_id] = mapped_reads
    return reads


def get_lineage(taxid, parent_map):
    lineage = []
    while taxid != "1" and taxid in parent_map:
        lineage.append(taxid)
        taxid = parent_map[taxid]
    lineage.append("1")
    return lineage


def summarize_krona(krona_file, parent_map, rank_map, diamond_reads=None):
    contig_counts = defaultdict(int)
    read_counts = defaultdict(int)
    totals_per_rank = defaultdict(int)

    with open(krona_file) as f:
        for line in f:
            if not line.strip():
                continue

            contig_id, taxid = line.strip().split("\t")
            lineage = get_lineage(taxid, parent_map)
            seen_ranks = set()

            for tid in lineage:
                r = rank_map.get(tid)
                if r in RANKS_OF_INTEREST and r not in seen_ranks:
                    key = (r, tid)
                    contig_counts[key] += 1
                    totals_per_rank[r] += 1

                    if diamond_reads is not None:
                        read_counts[key] += diamond_reads.get(contig_id, 0)

                    seen_ranks.add(r)

    return contig_counts, read_counts, totals_per_rank


def run(
    krona,
    diamond_tax,
    taxdump_nodes,
    taxdump_names,
    sample,
    classifier,
    unit,
    output,
):
    parent_map, rank_map, name_map = load_taxdump(taxdump_nodes, taxdump_names)

    diamond_reads = None
    if diamond_tax:
        diamond_reads = load_diamond_reads(diamond_tax)

    contig_counts, read_counts, totals_per_rank = summarize_krona(
        krona, parent_map, rank_map, diamond_reads
    )

    header = [
        "sample",
        "tool",
        "mode",
        "rank",
        "taxid",
        "name",
        "count",
        "percent",
        "source",
    ]

    if diamond_reads is not None:
        header.append("mapped_reads")

    with open(output, "w") as out:
        out.write("\t".join(header) + "\n")

        for rank in RANKS_OF_INTEREST:
            for (r, taxid), count in sorted(contig_counts.items()):
                if r != rank:
                    continue

                total = totals_per_rank[r]
                percent = (count / total * 100) if total else 0.0
                taxname = name_map.get(taxid, "NA")

                fields = [
                    sample,
                    classifier,
                    unit,
                    r,
                    taxid,
                    taxname,
                    str(count),
                    f"{percent:.4f}",
                    krona,
                ]

                if diamond_reads is not None:
                    fields.append(str(read_counts.get((r, taxid), 0)))

                out.write("\t".join(fields) + "\n")


# ---- Snakemake entrypoint ----
if "snakemake" in globals():
    run(
        krona=snakemake.input.krona,
        diamond_tax=snakemake.input.get("annotated", None),
        taxdump_nodes=os.path.join(snakemake.params.taxdump, "nodes.dmp"),
        taxdump_names=os.path.join(snakemake.params.taxdump, "names.dmp"),
        sample=snakemake.params.sample,
        classifier=snakemake.params.tool,
        unit=snakemake.params.mode,
        output=snakemake.output[0],
    )
else:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--krona", required=True)
    parser.add_argument("--diamond-tax")
    parser.add_argument("--taxdump-dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--classifier", required=True)
    parser.add_argument("--unit", default="contigs")
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    run(
        krona=args.krona,
        diamond_tax=args.diamond_tax,
        taxdump_nodes=os.path.join(args.taxdump_dir, "nodes.dmp"),
        taxdump_names=os.path.join(args.taxdump_dir, "names.dmp"),
        sample=args.sample,
        classifier=args.classifier,
        unit=args.unit,
        output=args.output,
    )
