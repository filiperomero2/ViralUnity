#!/usr/bin/env python

import os
from collections import defaultdict

def load_assembly_taxid_map(assembly_summary_path):
    assembly2taxid = {}
    with open(assembly_summary_path) as f:
        for line in f:
            if not line.startswith("#"):
                cols = line.strip().split("\t")
                if len(cols) > 5:
                    assembly_acc = cols[0].strip()
                    taxid = cols[5].strip()
                    assembly2taxid[assembly_acc] = taxid
    return assembly2taxid

def parse_diamond_results(diamond_output_path, assembly2taxid):
    read2taxid = {}
    # If diamond output is empty, return empty taxid mapping
    if os.stat(diamond_output_path).st_size == 0:
        return read2taxid
    with open(diamond_output_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            read_id = parts[0]
            sseqid = parts[1]
            assembly_acc = sseqid.split("|")[0]
            prot_acc = sseqid.split("|")[1]
            if assembly_acc == "":
                print(f"No assembly accession identified for {prot_acc}, treating as unclassified...")
            taxid = assembly2taxid.get(assembly_acc)
            if taxid:
                read2taxid[read_id] = taxid
    return read2taxid

def extract_ids_from_sequences(sequences,data_format):
    ids = set()
    if os.stat(sequences).st_size == 0:
        return ids
    with open(sequences) as f:
        if data_format == "fastq":
            for i, line in enumerate(f):
                if i % 4 == 0:  # fastq header
                    id = line.strip().lstrip("@").split()[0]
                    ids.add(id)
        else:
            for i, line in enumerate(f):
                if line.startswith(">"): # fasta header
                    id = line.strip().lstrip(">").split()[0]
                    ids.add(id)
    return ids

def write_krona_input(output_path, ids, read2taxid):
    with open(output_path, "w") as out:
        if len(ids) > 0:
            for id in ids:
                taxid = read2taxid.get(id, "0")
                out.write(f"{id}\t{taxid}\n")
        else:
            print(f"An empty krona input file was created -> {output_path}\n")

def main(input,output,data_format):
    diamond_output = input[0]
    sequences = input[1]
    assembly_summary = input[2]
    assembly2taxid = load_assembly_taxid_map(assembly_summary)
    read2taxid = parse_diamond_results(diamond_output, assembly2taxid)
    ids = extract_ids_from_sequences(sequences,data_format)
    write_krona_input(output, ids, read2taxid)
    print(f"Krona input written to {output}")

if __name__ == "__main__":
    input = snakemake.input
    output = snakemake.output[0]
    data_format = snakemake.params[0]
    main(input,output,data_format)
    exit()
