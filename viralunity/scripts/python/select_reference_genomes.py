#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import glob
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select reference genomes based on taxonomy or similarity."
    )
    parser.add_argument(
        "--summary-dir", required=True, help="Directory containing TSV summary files."
    )
    parser.add_argument(
        "--method", required=True, choices=["kraken2", "diamond", "both"]
    )
    parser.add_argument("--source", required=True, choices=["reads", "contigs", "both"])
    parser.add_argument("--reads-count", type=int, default=100)
    parser.add_argument("--contigs-count", type=int, default=1)
    parser.add_argument(
        "--families", required=True, help="Comma-separated list of families."
    )
    parser.add_argument("--strategy", required=True, choices=["taxid", "similarity"])
    parser.add_argument(
        "--genome2taxid", required=True, help="Path to genome2taxid.tsv"
    )
    parser.add_argument(
        "--blast-db", required=True, help="Path to viral.genomes.fasta BLAST DB"
    )
    parser.add_argument("--blast-qcov", type=float, default=80)
    parser.add_argument("--blast-pident", type=float, default=80)
    parser.add_argument(
        "--contigs-dir",
        help="Path to denovo assembly contigs directory (for similarity strategy)",
    )
    parser.add_argument("--out-tsv", required=True, help="Output TSV mapping file")
    return parser.parse_args()


def main(args=None):
    if args is None:
        args = parse_args()
    families = [f.strip() for f in args.families.split(",")]

    # Parse summary TSVs
    summary_files = []
    if args.method in ["kraken2", "both"]:
        if args.source in ["reads", "both"]:
            summary_files.extend(
                glob.glob(
                    os.path.join(args.summary_dir, "kraken2_reads", "kraken2_reads_taxa_summary.tsv")
                )
            )
        if args.source in ["contigs", "both"]:
            summary_files.extend(
                glob.glob(
                    os.path.join(args.summary_dir, "kraken2_contigs", "kraken2_contigs_taxa_summary.tsv")
                )
            )

    if args.method in ["diamond", "both"]:
        if args.source in ["reads", "both"]:
            summary_files.extend(
                glob.glob(
                    os.path.join(args.summary_dir, "diamond_reads", "diamond_reads_taxa_summary.tsv")
                )
            )
        if args.source in ["contigs", "both"]:
            summary_files.extend(
                glob.glob(
                    os.path.join(args.summary_dir, "diamond_contigs", "diamond_contigs_taxa_summary.tsv")
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

    # Filter by count
    df_reads = df[df["mode"] == "reads"]
    df_contigs = df[df["mode"] == "contigs"]

    df_reads = df_reads[df_reads["count"] >= args.reads_count]
    df_contigs = df_contigs[df_contigs["count"] >= args.contigs_count]

    filtered_df = pd.concat([df_reads, df_contigs])

    # Check which samples have the requested families
    family_hits = filtered_df[filtered_df["name"].isin(families)]
    samples_families = family_hits.groupby("sample")["name"].unique().to_dict()

    valid_samples = list(samples_families.keys())

    out_records = []

    if args.strategy == "taxid":
        # Load genome2taxid
        if not os.path.exists(args.genome2taxid):
            print(f"Warning: {args.genome2taxid} not found!")
            sys.exit(0)

        g2t = pd.read_csv(
            args.genome2taxid, sep="\t", header=None, names=["accession", "taxid"]
        )
        g2t["taxid"] = g2t["taxid"].astype(str)

        # Use all ranks (species/genus/family) for valid samples so species-level
        # taxids can match against genome2taxid (which holds species/strain taxids)
        taxid_candidates = filtered_df[filtered_df["sample"].isin(valid_samples)]
        taxid_candidates["taxid_str"] = taxid_candidates["taxid"].astype(str)

        merged = pd.merge(
            taxid_candidates, g2t, left_on="taxid_str", right_on="taxid", how="inner"
        )

        for idx, row in merged.iterrows():
            sample = row["sample"]
            acc = row["accession"]
            # ref_key uniquely identifies this assembly target: family_accession (sorted for determinism)
            fams = "_".join(sorted(samples_families[sample])).replace(" ", "_")
            ref_key = f"{fams}_{acc}"
            out_records.append(
                {"sample": sample, "ref_key": ref_key, "reference_genome": acc}
            )

    elif args.strategy == "similarity":
        # Blast contigs of valid samples against viral.genomes.fasta
        if not args.contigs_dir:
            print("Error: --contigs-dir required for similarity strategy.")
            sys.exit(1)

        for sample in valid_samples:
            contig_file = os.path.join(args.contigs_dir, f"sample-{sample}.viral_contigs.fa")
            if not os.path.exists(contig_file):
                contig_file_alt = os.path.join(args.contigs_dir, f"{sample}.viral_contigs.fa")
                if os.path.exists(contig_file_alt):
                    contig_file = contig_file_alt
                else:
                    contig_globs = glob.glob(os.path.join(args.contigs_dir, sample, "*.contigs.fa*"))
                    if not contig_globs:
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

            try:
                res = subprocess.run(blast_cmd, capture_output=True, text=True)
                if res.returncode == 0 and res.stdout:
                    for line in res.stdout.strip().split("\\n"):
                        parts = line.split("\\t")
                        if len(parts) >= 4:
                            sseqid = parts[1]
                            pident = float(parts[2])
                            qcov = float(parts[3])
                            if pident >= args.blast_pident and qcov >= args.blast_qcov:
                                fams = "_".join(sorted(samples_families[sample])).replace(
                                    " ", "_"
                                )
                                ref_key = f"{fams}_{sseqid}"
                                out_records.append(
                                    {
                                        "sample": sample,
                                        "ref_key": ref_key,
                                        "reference_genome": sseqid,
                                    }
                                )
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
    args.out_tsv = snakemake.output.tsv
    
    main(args)
elif __name__ == "__main__":
    main()
