# Metagenomics Nanopore v2 Workflow

This workflow performs **viral metagenomic analysis of Nanopore sequencing data**
using a contig-based strategy. It extends the original ViralUnity metagenomics
pipeline by adding:

1. Optional **host read removal (dehosting)**
2. **De novo assembly** of reads into contigs
3. **Protein-level taxonomic classification** using DIAMOND
4. Optional **polishing of viral contigs** using Medaka
5. Cohort-level **taxonomic summaries and Krona visualizations**

The workflow is implemented in Snakemake and is designed to be executed
either directly with Snakemake or via the ViralUnity CLI.

---

## High-level workflow overview

```
Nanopore reads
   |
   |-- (optional) Host read removal (minimap2 + samtools)
   |
   |-- De novo assembly (MEGAHIT)
   |
   |-- Taxonomic classification (both diamond and kraken2)
   |
   |-- Trimming of viral contig regions (DIAMOND-supported)
   |
   |-- Polishing of viral contigs (Medaka)
   |
   |-- Read-support filtering (samtools idxstats)
   |
   |-- Krona visualizations
   |
   `-- Summary tables
```

---

## Inputs

### Required inputs

- Nanopore sequencing reads (FASTQ)
- Sample sheet (provided via config)
- DIAMOND protein database (FASTA; indexed automatically)
- Krona taxonomy database

### Optional inputs

- **Host reference genome (FASTA)**  
  If provided, host-associated reads are removed prior to assembly.

- **Kraken2 database**  
  Used only if `run_kraken2: true`.

- **Taxonomy resources** (required for DIAMOND annotation):
  - NCBI taxdump
  - Assembly summary file
  - TaxID-to-family mapping file

---

## Configuration

A minimal example configuration file:

```yaml
samples:
  sample-barcode05: /path/to/barcode05.fastq
  sample-barcode09: /path/to/barcode09.fastq

data: nanopore
output: output/meta_nanopore_run/
threads: 2

# Pipeline toggles
run_denovo_assembly: true
run_kraken2: false
run_diamond: true

# Host filtering
host_reference: "NA"   # or path to host FASTA

# Databases
kraken2_database: "/path/to/kraken2_db/"
krona_database: "/path/to/krona/taxonomy/"
diamond_database: "/path/to/diamond_db.fasta"

diamond_sensitivity: "sensitive"
evalue: 1e-10

# Taxonomy resources
taxdump: "/path/to/taxdump/"
assembly_summary: "/path/to/assembly_summary.tsv"
taxid_to_family: "/path/to/taxid_to_family.csv"

# Filtering options
remove_human_reads: false
remove_unclassified_reads: false

# Medaka
medaka_model: "r941_min_high_g360"
```

A template configuration is provided in `config.default.yml`.

---

## Running the workflow

### Using Snakemake directly

From the workflow directory:

```bash
snakemake   --snakefile Snakefile   --configfile config.test.yml   --cores 4   -p
```
---

## Outputs

All outputs are written under the configured `output/` directory.

### De novo assembly
```
denovo_assembly/megahit/{sample}/final.contigs.fa
```

### DIAMOND contig classification
```
metagenomics/taxonomic_assignments/diamond_contigs/results/
  |-- {sample}.diamond.tsv
  |-- {sample}.diamond.supported.tsv
  |-- {sample}.diamond.supported.tax.tsv
  `-- {sample}.diamond.trimmed.fa
```

### Polished viral contigs (Medaka)
```
medaka/{sample}/viral_consensus.fasta
```

### Krona reports
```
metagenomics/taxonomic_assignments/diamond_contigs/reports/
  `-- {sample}.diamond.supported.krona.html
```

### Cohort-level summaries
```
metagenomics/taxonomic_assignments/diamond_contigs/
  `-- diamond_contigs_metagenomics_summary.txt
```

### Logs and benchmarks
```
logs/
  |-- megahit/
  |-- diamond_contigs/
  |-- medaka_consensus_trimmed/
  `-- ...
```

---

## Notes and design considerations

- DIAMOND classification is performed on **assembled contigs**, not raw reads.
- Polishing with Medaka is applied only to **DIAMOND-supported viral regions**.
- Krona visualizations are generated only when valid input data is available.

---

## Status

This workflow represents an **extended metagenomics pipeline** for Nanopore
data in ViralUnity. At this time, it is intended to coexist with the original metagenomics
workflow and may evolve further as additional classifiers and post-processing
steps are added.
