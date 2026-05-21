#!/usr/bin/env python3
import argparse
import glob
import os
import re
import subprocess
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select reference genomes based on taxonomy or similarity."
    )
    parser.add_argument(
        "--summary-dir", required=True, help="Directory containing TSV summary files."
    )
    parser.add_argument("--method", required=True, choices=["kraken2", "diamond", "both"])
    parser.add_argument("--source", required=True, choices=["reads", "contigs", "both"])
    parser.add_argument("--reads-count", type=int, default=100)
    parser.add_argument("--contigs-count", type=int, default=1)
    parser.add_argument("--families", required=True, help="Comma-separated list of families.")
    parser.add_argument("--strategy", required=True, choices=["taxid", "similarity"])
    parser.add_argument("--genome2taxid", required=True, help="Path to genome2taxid.tsv")
    parser.add_argument("--blast-db", required=True, help="Path to viral.genomes.fasta BLAST DB")
    parser.add_argument("--blast-qcov", type=float, default=80)
    parser.add_argument("--blast-pident", type=float, default=80)
    parser.add_argument(
        "--contigs-dir",
        help="Path to denovo assembly contigs directory (for similarity strategy)",
    )
    parser.add_argument(
        "--taxdump",
        default="",
        help="Path to NCBI taxdump directory (nodes.dmp + names.dmp) used to resolve "
        "the viral family for each matched accession.",
    )
    parser.add_argument("--out-tsv", required=True, help="Output TSV mapping file")
    return parser.parse_args()


def load_taxdump(taxdump_dir):
    """Return (nodes, names) dicts loaded from NCBI taxdump files.

    nodes : taxid_str -> (parent_taxid_str, rank_str)
    names : taxid_str -> scientific_name_str
    """
    nodes = {}
    names = {}
    nodes_file = os.path.join(taxdump_dir, "nodes.dmp")
    names_file = os.path.join(taxdump_dir, "names.dmp")
    with open(nodes_file) as f:
        for line in f:
            parts = line.split("|")
            taxid = parts[0].strip()
            parent = parts[1].strip()
            rank = parts[2].strip()
            nodes[taxid] = (parent, rank)
    with open(names_file) as f:
        for line in f:
            parts = line.split("|")
            taxid = parts[0].strip()
            name = parts[1].strip()
            name_class = parts[3].strip()
            if name_class == "scientific name":
                names[taxid] = name
    return nodes, names


def get_family_for_taxid(taxid_str, nodes, names):
    """Walk the taxonomy tree upward and return the family-rank name, or None."""
    visited = set()
    current = taxid_str
    while current and current not in visited:
        visited.add(current)
        if current not in nodes:
            break
        parent, rank = nodes[current]
        if rank == "family":
            return names.get(current)
        if current == parent:  # reached taxonomy root
            break
        current = parent
    return None


def get_taxid_at_rank(taxid_str, target_rank, nodes):
    """Walk the taxonomy tree upward and return the taxid string at target_rank, or None."""
    visited = set()
    current = taxid_str
    while current and current not in visited:
        visited.add(current)
        if current not in nodes:
            break
        parent, rank = nodes[current]
        if rank == target_rank:
            return current
        if current == parent:  # reached taxonomy root
            break
        current = parent
    return None


def clean_accession(sseqid):
    """Strip NCBI pipe-delimited prefixes (e.g. ref|NC_004911.1|) from a BLAST sseqid."""
    m = re.match(r"^[a-z]+\|([^|]+)\|?$", sseqid)
    if m:
        return m.group(1)
    return sseqid


def main(args=None):
    if args is None:
        args = parse_args()
    families = [f.strip() for f in args.families.split(",")]

    # Both strategies require a populated genome2taxid mapping: taxid uses it
    # as the primary accession lookup; similarity uses it to validate that
    # BLAST hits trace back to a target family. A missing or empty file
    # would otherwise yield a silent header-only reference_targets.tsv.
    if not os.path.exists(args.genome2taxid):
        print(
            f"ERROR: --viral-taxids path does not exist: {args.genome2taxid}\n"
            "Run 'viralunity get-databases virus-genome' to download the viral "
            "genomes database and build the taxid mapping.",
            file=sys.stderr,
        )
        sys.exit(1)
    if os.path.getsize(args.genome2taxid) == 0:
        print(
            f"ERROR: --viral-taxids file is empty: {args.genome2taxid}\n"
            "Re-run 'viralunity get-databases virus-genome' to (re)populate it.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Load taxdump for family-aware filtering (optional but strongly recommended)
    nodes, tax_names = {}, {}
    taxdump_dir = getattr(args, "taxdump", "") or ""
    if taxdump_dir and os.path.isdir(taxdump_dir):
        print(f"Loading taxdump from: {taxdump_dir}")
        nodes, tax_names = load_taxdump(taxdump_dir)
        print(f"  Loaded {len(nodes)} taxonomy nodes.")
    else:
        print(
            "Warning: --taxdump not provided or directory not found. "
            "Family-level filtering will be skipped and ref_key may be inaccurate."
        )

    # Parse summary TSVs
    summary_files = []
    if args.method in ["kraken2", "both"]:
        if args.source in ["reads", "both"]:
            summary_files.extend(
                glob.glob(
                    os.path.join(
                        args.summary_dir, "kraken2_reads", "kraken2_reads_taxa_summary.tsv"
                    )
                )
            )
        if args.source in ["contigs", "both"]:
            summary_files.extend(
                glob.glob(
                    os.path.join(
                        args.summary_dir, "kraken2_contigs", "kraken2_contigs_taxa_summary.tsv"
                    )
                )
            )

    if args.method in ["diamond", "both"]:
        if args.source in ["reads", "both"]:
            summary_files.extend(
                glob.glob(
                    os.path.join(
                        args.summary_dir, "diamond_reads", "diamond_reads_taxa_summary.tsv"
                    )
                )
            )
        if args.source in ["contigs", "both"]:
            summary_files.extend(
                glob.glob(
                    os.path.join(
                        args.summary_dir, "diamond_contigs", "diamond_contigs_taxa_summary.tsv"
                    )
                )
            )

    if not summary_files:
        print("Warning: No summary TSVs found for the given method and source.")
        pd.DataFrame(columns=["sample", "ref_key", "reference_genome"]).to_csv(
            args.out_tsv, sep="\t", index=False
        )
        return

    frames = []
    for f in summary_files:
        if os.path.exists(f):
            frames.append(pd.read_csv(f, sep="\t"))

    if not frames:
        pd.DataFrame(columns=["sample", "ref_key", "reference_genome"]).to_csv(
            args.out_tsv, sep="\t", index=False
        )
        return

    df = pd.concat(frames, ignore_index=True)

    # Filter by count threshold
    df_reads = df[df["mode"] == "reads"]
    df_contigs = df[df["mode"] == "contigs"]

    df_reads = df_reads[df_reads["count"] >= args.reads_count]
    df_contigs = df_contigs[df_contigs["count"] >= args.contigs_count]

    filtered_df = pd.concat([df_reads, df_contigs])

    # Identify samples that have hits to any of the requested families
    family_hits = filtered_df[filtered_df["name"].isin(families)]
    valid_samples = list(family_hits["sample"].unique())

    out_records = []

    if args.strategy == "taxid":
        g2t = pd.read_csv(args.genome2taxid, sep="\t", header=None, names=["accession", "taxid"])
        g2t["taxid"] = g2t["taxid"].astype(str)
        g2t["accession"] = g2t["accession"].astype(str)

        # Build a direct taxid -> accessions map.
        # The taxid value from the summary is used as-is (no rank normalisation).
        # This avoids cross-database taxid mismatches that silently produce empty results
        # when protein2taxid and genome2taxid use different taxid levels for the same organism.
        taxid_to_accs = {}
        for _, grow in g2t.iterrows():
            taxid_to_accs.setdefault(grow["taxid"], []).append(grow["accession"])

        taxid_candidates = filtered_df[filtered_df["sample"].isin(valid_samples)].copy()
        taxid_candidates["taxid_str"] = taxid_candidates["taxid"].astype(str)

        # Track which (sample, family) pairs produced at least one accession so we can
        # emit a single warning per family rather than one noisy warning per taxid row.
        found_per_sample_family = {}

        for _, row in taxid_candidates.iterrows():
            sample = row["sample"]
            taxid_str = row["taxid_str"]

            # Use taxdump only to validate the hit belongs to a target family and to
            # build a clean ref_key — not to normalise the taxid lookup level.
            if nodes:
                family = get_family_for_taxid(taxid_str, nodes, tax_names)
                if family is None or family not in families:
                    continue
            else:
                all_fams = family_hits[family_hits["sample"] == sample]["name"].unique()
                family = sorted(all_fams)[0] if len(all_fams) else "unknown"

            accs = taxid_to_accs.get(taxid_str, [])
            if not accs and nodes:
                # Fallback: if the exact taxid isn't in genome2taxid (e.g. it's a
                # strain or subspecies), try the species-level ancestor taxid.
                species_taxid = get_taxid_at_rank(taxid_str, "species", nodes)
                if species_taxid and species_taxid != taxid_str:
                    accs = taxid_to_accs.get(species_taxid, [])
            if not accs:
                continue  # silently skip; will warn per-family below if nothing found

            found_per_sample_family.setdefault(sample, set()).add(family)
            for acc in accs:
                ref_key = f"{family.replace(' ', '_')}_{acc}"
                out_records.append({"sample": sample, "ref_key": ref_key, "reference_genome": acc})

        # Warn once per (sample, target-family) that produced no accessions at all.
        for sample in valid_samples:
            sample_families = set(
                family_hits[family_hits["sample"] == sample]["name"].unique()
            ) & set(families)
            for fam in sample_families:
                if fam not in found_per_sample_family.get(sample, set()):
                    print(
                        f"  Warning: no genome found in genome2taxid for family {fam!r}, "
                        f"sample {sample}. If using a --refseq database, only one accession "
                        "per species is indexed; for broader coverage use "
                        "--reference-selection-strategy similarity."
                    )

    elif args.strategy == "similarity":
        if not args.contigs_dir:
            print("Error: --contigs-dir required for similarity strategy.")
            sys.exit(1)

        # Load genome2taxid for family lookup of BLAST hits (taxdump not required)
        g2t = pd.read_csv(
            args.genome2taxid, sep="\t", header=None, names=["accession", "taxid"]
        )
        g2t_dict = dict(zip(g2t["accession"].astype(str), g2t["taxid"].astype(str)))

        for sample in valid_samples:
            contig_file = os.path.join(args.contigs_dir, f"sample-{sample}.viral_contigs.fa")
            if not os.path.exists(contig_file):
                contig_file_alt = os.path.join(args.contigs_dir, f"{sample}.viral_contigs.fa")
                if os.path.exists(contig_file_alt):
                    contig_file = contig_file_alt
                else:
                    contig_globs = glob.glob(
                        os.path.join(args.contigs_dir, sample, "*.contigs.fa*")
                    )
                    if not contig_globs:
                        print(
                            f"Warning: no contig file found for sample {sample} in {args.contigs_dir}, skipping."
                        )
                        continue
                    contig_file = contig_globs[0]

            blast_cmd = [
                "blastn",
                "-query",
                contig_file,
                "-db",
                args.blast_db,
                "-outfmt",
                "6 qseqid sseqid pident qcovhsp",
                "-max_target_seqs",
                "10",
            ]

            print(f"Running BLAST for {sample}: {' '.join(blast_cmd)}")
            try:
                res = subprocess.run(blast_cmd, capture_output=True, text=True)
                if res.returncode != 0:
                    print(f"BLAST failed for {sample} (exit {res.returncode}):\n{res.stderr}")
                    continue
                if not res.stdout.strip():
                    print(f"BLAST produced no hits for {sample}.")
                    continue

                # BLAST output is sorted by bitscore (best first) within each query.
                # Track assigned contigs so only the best qualifying hit per contig is kept.
                assigned_contigs = set()
                for line in res.stdout.strip().split("\n"):
                    parts = line.split("\t")
                    if len(parts) < 4:
                        continue
                    qseqid = parts[0]
                    if qseqid in assigned_contigs:
                        continue
                    raw_sseqid = parts[1]
                    try:
                        pident = float(parts[2])
                        qcov = float(parts[3])
                    except ValueError:
                        continue
                    if pident < args.blast_pident or qcov < args.blast_qcov:
                        continue

                    # Strip NCBI pipe-delimited prefix (e.g. ref|NC_004911.1|)
                    acc = clean_accession(raw_sseqid)

                    # Keep hits whose taxid (at any rank) traces back to a target family.
                    if g2t_dict:
                        hit_taxid = g2t_dict.get(acc, "")
                        if not hit_taxid:
                            continue
                        if nodes:
                            family = get_family_for_taxid(hit_taxid, nodes, tax_names)
                            if family is None or family not in families:
                                continue
                            ref_key = f"{family.replace(' ', '_')}_{acc}"
                        else:
                            all_fams = family_hits[family_hits["sample"] == sample]["name"].unique()
                            fams = "_".join(sorted(all_fams)).replace(" ", "_")
                            ref_key = f"{fams}_{acc}"
                    else:
                        all_fams = family_hits[family_hits["sample"] == sample]["name"].unique()
                        fams = "_".join(sorted(all_fams)).replace(" ", "_")
                        ref_key = f"{fams}_{acc}"

                    out_records.append(
                        {"sample": sample, "ref_key": ref_key, "reference_genome": acc}
                    )
                    assigned_contigs.add(qseqid)
            except Exception as e:
                print(f"Error running BLAST for {sample}: {e}")

    # Deduplicate and save
    if out_records:
        out_df = pd.DataFrame(out_records).drop_duplicates()
    else:
        out_df = pd.DataFrame(columns=["sample", "ref_key", "reference_genome"])

    out_df.to_csv(args.out_tsv, sep="\t", index=False)


if "snakemake" in globals():

    class SnakemakeArgs:
        pass

    args = SnakemakeArgs()
    args.summary_dir = snakemake.params.summary_dir
    args.method = snakemake.params.method
    args.source = snakemake.params.source
    args.reads_count = snakemake.params.reads_count
    args.contigs_count = snakemake.params.contigs_count
    args.families = snakemake.params.families
    args.strategy = snakemake.params.strategy
    args.genome2taxid = snakemake.params.genome2taxid
    args.blast_db = snakemake.params.blast_db
    args.blast_qcov = snakemake.params.blast_qcov
    args.blast_pident = snakemake.params.blast_pident
    args.contigs_dir = snakemake.params.contigs_dir
    args.taxdump = snakemake.params.taxdump
    args.out_tsv = snakemake.output.tsv

    main(args)
elif __name__ == "__main__":
    main()
