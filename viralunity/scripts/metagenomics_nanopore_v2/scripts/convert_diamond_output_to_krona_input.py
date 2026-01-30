#!/usr/bin/env python3

import os
import gzip


def open_maybe_gzip(path, mode="rt"):
    """
    Open plain-text or gzipped text file.
    mode should be 'rt' or 'wt' for text, 'rb'/'wb' for binary.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def load_assembly_taxid_map(assembly_summary_path):
    assembly2taxid = {}
    with open_maybe_gzip(assembly_summary_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 5:
                assembly_acc = cols[0].strip()
                taxid = cols[5].strip()
                if assembly_acc and taxid:
                    assembly2taxid[assembly_acc] = taxid
    return assembly2taxid


def parse_diamond_results(diamond_output_path, assembly2taxid):
    read2taxid = {}

    # If file is truly empty on disk, skip quickly
    try:
        if os.path.getsize(diamond_output_path) == 0:
            return read2taxid
    except OSError:
        return read2taxid

    with open_maybe_gzip(diamond_output_path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            read_id = parts[0]
            sseqid = parts[1]

            # Expected: assembly|protein (your previous format)
            fields = sseqid.split("|")
            assembly_acc = fields[0] if len(fields) >= 1 else ""
            prot_acc = fields[1] if len(fields) >= 2 else ""

            if not assembly_acc:
                # Keep unclassified by omission; krona will get taxid 0 later
                # (avoid noisy prints in big runs)
                continue

            taxid = assembly2taxid.get(assembly_acc)
            if taxid:
                read2taxid[read_id] = taxid

    return read2taxid


def extract_ids_from_sequences(sequences, data_format):
    ids = set()

    # Avoid relying on compressed "emptiness"; just iterate safely.
    with open_maybe_gzip(sequences, "rt") as f:
        if data_format == "fastq":
            for i, line in enumerate(f):
                if i % 4 == 0:  # fastq header
                    rid = line.strip().lstrip("@").split()[0]
                    if rid:
                        ids.add(rid)
        else:
            for line in f:
                if line.startswith(">"):  # fasta header
                    rid = line.strip().lstrip(">").split()[0]
                    if rid:
                        ids.add(rid)

    return ids


def write_krona_input(output_path, ids, read2taxid):
    with open(output_path, "w") as out:
        if ids:
            for rid in ids:
                taxid = read2taxid.get(rid, "0")
                out.write(f"{rid}\t{taxid}\n")
        else:
            # keep file empty but valid
            pass


def main(input_files, output_file, data_format):
    diamond_output = input_files[0]
    sequences = input_files[1]
    assembly_summary = input_files[2]

    assembly2taxid = load_assembly_taxid_map(assembly_summary)
    read2taxid = parse_diamond_results(diamond_output, assembly2taxid)
    ids = extract_ids_from_sequences(sequences, data_format)
    write_krona_input(output_file, ids, read2taxid)

    print(f"Krona input written to {output_file}")


if __name__ == "__main__":
    # snakemake object is injected by Snakemake when using `script:`
    input_files = list(snakemake.input)
    output_file = snakemake.output[0]

    # Prefer named param if present; fallback to positional
    try:
        data_format = snakemake.params.data_format
    except Exception:
        data_format = snakemake.params[0]

    main(input_files, output_file, data_format)
