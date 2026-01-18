#!/usr/bin/env python3
"""Annotate a DIAMOND outfmt-6 file (already "supported" & read‑count annotated)
with taxonomic family / genus / species using:

* NCBI **viral RefSeq assembly summary** (provides assembly accession → taxid)
* NCBI **taxdump** (`nodes.dmp`, `names.dmp`) to climb the lineage

The script **does not** require internet and is fast once the look‑up tables are
in memory.

───────────────────────────────────────────────────────────────────────────────
Expected columns in the input DIAMOND file (minimum):
    qseqid   sseqid   ...   mapped_reads
The script only reads **sseqid**, appends 3 columns:
    family   genus   species

`sseqid` format example (viral RefSeq):
    GCF_000848265.1|NP_955595.1|ABL|Abelson_murine_leukemia_virus
The first pipe‑delimited field is the assembly accession.

───────────────────────────────────────────────────────────────────────────────
CLI example
~~~~~~~~~~~
```
python annotate_diamond_taxonomy.py \
  --diamond  sample.diamond.supported.tsv \
  --assembly viral_refseq_assembly_summary.tsv \
  --taxdump  taxdump/ \
  --output   sample.diamond.supported.tax.tsv
```

Snakemake rule snippet
~~~~~~~~~~~~~~~~~~~~~~
```
rule annotate_diamond_taxonomy:
    input:
        diamond   = "results/{sample}.diamond.supported.tsv",
        assembly  = "resources/viral_refseq_assembly_summary.tsv",
        taxdump   = directory("resources/taxdump")  # contains nodes.dmp, names.dmp
    output:
        annotated = "results/{sample}.diamond.supported.tax.tsv"
    log:
        "logs/annotate_tax/{sample}.log"
    script:
        "workflow/scripts/annotate_diamond_taxonomy.py"
```
"""

import argparse
import csv
import gzip
from pathlib import Path
from typing import Dict, Tuple, Optional

###############################################################################
# Argument parsing                                                             
###############################################################################

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Annotate DIAMOND hits with taxonomy")
    p.add_argument("--diamond", required=True, help="Input supported DIAMOND TSV")
    p.add_argument("--assembly", required=True, help="viral_refseq_assembly_summary.tsv")
    p.add_argument("--taxdump", required=True, help="Directory with nodes.dmp & names.dmp")
    p.add_argument("--output", required=True, help="Output annotated TSV path")
    return p.parse_args()

###############################################################################
# Taxdump parsing                                                              
###############################################################################

_TaxName = Tuple[str, str]  # (scientific_name, rank)


def load_names(names_path: Path) -> Dict[int, str]:
    """Return {taxid: scientific_name}. Only keeps names with class 'scientific name'."""
    names: Dict[int, str] = {}
    opener = gzip.open if names_path.suffix == ".gz" else open
    with opener(names_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            taxid, name_txt, _, name_class = parts[:4]
            if name_class == "scientific name":
                try:
                    names[int(taxid)] = name_txt
                except ValueError:
                    continue
    return names


def load_nodes(nodes_path: Path) -> Dict[int, Tuple[int, str]]:
    """Return {taxid: (parent_taxid, rank)}."""
    nodes: Dict[int, Tuple[int, str]] = {}
    opener = gzip.open if nodes_path.suffix == ".gz" else open
    with opener(nodes_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            taxid, parent, rank = parts[:3]
            try:
                nodes[int(taxid)] = (int(parent), rank)
            except ValueError:
                continue
    return nodes


def lineage_for(taxid: int, nodes: Dict[int, Tuple[int, str]], names: Dict[int, str]) -> Dict[str, str]:
    """Return mapping {rank: scientific_name} while climbing to root."""
    lineage: Dict[str, str] = {}
    current = taxid
    seen = set()
    while current not in seen and current in nodes:
        seen.add(current)
        parent, rank = nodes[current]
        if rank in {"species", "genus", "family"}:
            lineage[rank] = names.get(current, "NA")
            if len(lineage) == 3:
                break  # got all ranks
        if current == parent:
            break
        current = parent
    return lineage

###############################################################################
# Assembly summary parsing                                                     
###############################################################################

def load_assembly_taxids(assembly_path: Path) -> Dict[str, int]:
    """Return {assembly_accession (no version): taxid}.

    Handles the standard RefSeq assembly summary where the first header field is
    "# assembly_accession" (note the leading "# ").  Works for either un‑gzipped
    or gzipped (`.gz`) files.
    """
    mapping: Dict[str, int] = {}
    opener = gzip.open if assembly_path.suffix == ".gz" else open
    with opener(assembly_path, "rt", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="	")
        header = next(reader)
        # Clean the first header field (remove leading "# ") if present
        header[0] = header[0].lstrip("# ")
        try:
            idx_acc = header.index("assembly_accession")
            idx_tax = header.index("species_taxid")
        except ValueError as e:
            raise SystemExit("[ERROR] Could not find 'assembly_accession' or 'taxid' columns in assembly summary") from e

        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            acc = row[idx_acc].split(".")[0]
            try:
                mapping[acc] = int(row[idx_tax])
            except ValueError:
                continue  # skip malformed taxid
    return mapping

###############################################################################
# Annotation logic                                                             
###############################################################################

def annotate_diamond(
    diamond_path: Path,
    output_path: Path,
    acc_to_taxid: Dict[str, int],
    nodes: Dict[int, Tuple[int, str]],
    names: Dict[int, str],
) -> Tuple[int, int]:
    """Read DIAMOND file, write annotated version. Returns (rows_written, missing)."""
    written = 0
    missing = 0
    with diamond_path.open() as inp, output_path.open("w") as out:
        for line in inp:
            if line.startswith("#"):
                if line.startswith("# added_column: mapped_reads"):
                    out.write(line)  # keep existing annotation headers
                continue  # skip other diamond comments to avoid bloat
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            sseqid = parts[1]
            assembly_acc = sseqid.split("|", 1)[0].split(".")[0]
            taxid = acc_to_taxid.get(assembly_acc)
            family = genus = species = "NA"
            if taxid:
                lineage = lineage_for(taxid, nodes, names)
                family = lineage.get("family", "NA")
                genus = lineage.get("genus", "NA")
                species = lineage.get("species", "NA")
            else:
                missing += 1
            out.write("\t".join(parts + [family, genus, species]) + "\n")
            written += 1
    return written, missing

###############################################################################
# main()                                                                       
###############################################################################

def main() -> None:
    if "snakemake" in globals():
        diamond_path = Path(snakemake.input.diamond)
        assembly_path = Path(snakemake.input.assembly)
        taxdump_dir = Path(snakemake.params.taxdump)
        output_path = Path(snakemake.output.annotated)
    else:
        args = parse_args()
        diamond_path = Path(args.diamond)
        assembly_path = Path(args.assembly)
        taxdump_dir = Path(args.taxdump)
        output_path = Path(args.output)

    nodes_path = taxdump_dir / "nodes.dmp"
    names_path = taxdump_dir / "names.dmp"
    if not nodes_path.exists() or not names_path.exists():
        raise SystemExit("taxdump directory must contain nodes.dmp and names.dmp")

    print("[+] Loading taxonomy dumps…", flush=True)
    names = load_names(names_path)
    nodes = load_nodes(nodes_path)

    print("[+] Loading assembly summary…", flush=True)
    acc_to_taxid = load_assembly_taxids(assembly_path)

    print("[+] Annotating DIAMOND table…", flush=True)
    written, missing = annotate_diamond(diamond_path, output_path, acc_to_taxid, nodes, names)

    print(
        f"[✓] Wrote {written} rows to {output_path}.\n"
        f"    Taxonomy missing for {missing} rows (assembly accession not found).")


if __name__ == "__main__":
    main()
