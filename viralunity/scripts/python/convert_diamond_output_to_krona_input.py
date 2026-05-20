#!/usr/bin/env python3

import gzip
import os


def open_maybe_gzip(path, mode="rt"):
    """
    Open plain-text or gzipped text file.
    mode should be 'rt' or 'wt' for text, 'rb'/'wb' for binary.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def load_assembly_taxid_map(taxid_map_path):
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
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            read_id = parts[0]
            sseqid = parts[1]

            # Expected: assembly|protein (e.g. NCBI RefSeq format).
            # The protein component is not used downstream; only the
            # assembly accession is looked up for taxid mapping.
            fields = sseqid.split("|")
            assembly_acc = fields[0] if len(fields) >= 1 else ""

            if not assembly_acc:
                continue

            taxid = assembly2taxid.get(assembly_acc)
            if taxid:
                read2taxid[read_id] = taxid

    return read2taxid


def extract_ids_from_sequences(sequences, data_format):
    ids = set()

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


def main(input_files, output_file, data_format):
    diamond_output = input_files[0]
    sequences = input_files[1]
    taxids = input_files[2]

    assembly2taxid = load_assembly_taxid_map(taxids)
    read2taxid = parse_diamond_results(diamond_output, assembly2taxid)
    ids = extract_ids_from_sequences(sequences, data_format)
    write_krona_input(output_file, ids, read2taxid)

    print(f"Krona input written to {output_file}")


if __name__ == "__main__":
    if "snakemake" in globals():
        # Named input: diamond, fastq or fasta, assembly
        inp = snakemake.input
        seq_file = getattr(inp, "fastq", None) or getattr(inp, "fasta", None)
        input_files = [inp.diamond, seq_file, inp.assembly]
        output_file = getattr(snakemake.output, "krona_input", snakemake.output[0])
        data_format = snakemake.params.data_format
        main(input_files, output_file, data_format)
    else:

        # CLI usage
        parser = argparse.ArgumentParser()
        parser.add_argument("--diamond", required=True)
        parser.add_argument("--sequences", required=True)
        parser.add_argument("--taxids", required=True)
        parser.add_argument("--output", required=True)
        parser.add_argument("--data-format", choices=["fastq", "fasta"], default="fastq")
        args = parser.parse_args()
        main(
            [args.diamond, args.sequences, args.taxids],
            args.output,
            args.data_format,
        )
