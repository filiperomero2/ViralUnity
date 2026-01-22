# Metagenomics Nanopore v2 Workflow

This workflow performs **viral metagenomic analysis of Nanopore sequencing data**
within the ViralUnity framework.

It is an extended and more flexible successor of the original metagenomics
pipeline, supporting **read-based and contig-based taxonomic classification**,
optional **host read removal**, and **protein-level annotation using DIAMOND**.

The workflow is implemented in **Snakemake** and can be executed either:
- directly with Snakemake, or
- via the **ViralUnity CLI** (`viralunity meta --pipeline v2`)

---

## Supported analysis modes

### Read-based analyses
- Kraken2 on reads
- DIAMOND on reads

### Contig-based analyses
- De novo assembly (MEGAHIT)
- Kraken2 on contigs
- DIAMOND blastx on contigs
- Medaka polishing

---

## Minimal configuration example

```yaml
samples:
  sample-barcode05: /path/to/barcode05.fastq
  sample-barcode09: /path/to/barcode09.fastq

data: nanopore
output: output/meta_nanopore_run/
threads: 2

run_denovo_assembly: true
run_kraken2_reads: true
run_diamond: true

kraken2_database: /path/to/kraken2_db/
krona_database: /path/to/krona/taxonomy/
diamond_database: /path/to/diamond_db.fasta

taxdump: /path/to/taxdump/
assembly_summary: /path/to/assembly_summary.tsv
taxid_to_family: /path/to/taxid_to_family.csv
```

---

## Running with ViralUnity CLI

```bash
viralunity meta \
  --data-type nanopore \
  --pipeline v2 \
  --sample-sheet samples.csv \
  --output output \
  --run-name meta_nanopore_v2 \
  --run-kraken2-reads \
  --run-diamond \
  --run-denovo-assembly \
  --kraken2-database /path/to/kraken2_db \
  --diamond-database /path/to/diamond_db.fasta \
  --taxdump /path/to/taxdump \
  --assembly-summary /path/to/assembly_summary.tsv \
  --taxid-to-family /path/to/taxid_to_family.csv \
  --krona-database /path/to/krona/taxonomy
```

---

## Outputs

- Kraken2 reports and summaries (reads and/or contigs)
- DIAMOND classification tables
- Krona interactive HTML reports
- Polished viral contigs

---

## Status

This workflow is the **recommended metagenomics pipeline for Nanopore data**
in ViralUnity and is under active development.
