#!/usr/bin/env python3
"""Produce a per-contig taxonomy table (qseqid, taxid, name, mapped_reads) from
DIAMOND supported TSV and assembly/taxdump, for use by summarize_krona_taxa.
"""

import gzip
import os
from pathlib import Path


def open_maybe_gzip(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def load_assembly_taxid(taxid_map_path):
    """Load genome accession -> TaxID mapping from protein2taxid.tsv.

    The file is produced by 'viralunity get-databases diamond' and has two
    tab-separated columns (no header):
        genome_accession  TAB  taxid

    Lines beginning with '#' are treated as comments and skipped.
    """
    acc2taxid = {}
    with open_maybe_gzip(taxid_map_path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= 2:
                acc = cols[0].strip()
                taxid = cols[1].strip()
                if acc and taxid:
                    acc2taxid[acc] = taxid
    return acc2taxid


def load_taxid_names(taxdump_dir):
    names_path = os.path.join(taxdump_dir, "names.dmp")
    taxid2name = {}
    with open(names_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                taxid2name[parts[0]] = parts[1]
    return taxid2name


def run(diamond_tsv, taxids, taxdump_dir, output_path):
    acc2taxid = load_assembly_taxid(taxids)
    taxid2name = load_taxid_names(taxdump_dir)

    with open(diamond_tsv) as inp, open(output_path, "w") as out:
        for line in inp:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            qseqid = parts[0]
            sseqid = parts[1]
            mapped_reads = parts[-1] if len(parts) > 1 else "0"
            acc = sseqid.split("|")[0] if "|" in sseqid else sseqid
            taxid = acc2taxid.get(acc, "0")
            name = taxid2name.get(taxid, "NA")
            out.write(f"{qseqid}\t{taxid}\t{name}\t{mapped_reads}\n")


if __name__ == "__main__":
    if "snakemake" in globals():
        run(
            diamond_tsv=snakemake.input.diamond,
            taxids=snakemake.input.assembly,
            taxdump_dir=snakemake.params.taxdump,
            output_path=snakemake.output.annotated,
        )
    else:
        import argparse

        p = argparse.ArgumentParser()
        p.add_argument("--diamond", required=True)
        p.add_argument("--assembly", required=True)
        p.add_argument("--taxdump", required=True)
        p.add_argument("--output", required=True)
        a = p.parse_args()
        run(a.diamond, a.assembly, a.taxdump, a.output)
