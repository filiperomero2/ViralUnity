# Usage Overview

ViralUnity provides five main subcommands:

```
viralunity [--version]
├── build-deacon-index   Build a Deacon minimizer index from a FASTA file
├── create-samplesheet   Generate a sample-sheet CSV from a run directory
├── get-databases        Download and set up reference databases
│   ├── kraken2              Download a Kraken2 pre-built index
│   ├── krona                Set up the Krona taxonomy database
│   ├── taxdump              Download the NCBI taxdump
│   ├── diamond              Download RefSeq viral proteins and build a Diamond DB
│   ├── clean-protein-fasta  Strip nucleotide records from a protein FASTA
│   ├── virus-genome         Download viral genomes and build a BLAST index
│   ├── host-genome          Download a host genome using NCBI Datasets
│   ├── deacon-index         Download a pre-built Deacon minimizer index
│   └── all                  Run every download/index step above in sequence
├── consensus
│   ├── illumina          Reference-based consensus assembly for Illumina data
│   └── nanopore          Reference-based consensus assembly for Nanopore data
└── meta
    ├── illumina          Metagenomics pipeline for Illumina data
    └── nanopore          Metagenomics pipeline for Nanopore data
```

Use `--help` at any level for the full option list:

```bash
viralunity --help
viralunity consensus --help
viralunity consensus illumina --help
viralunity meta nanopore --help
```

## General Workflow

1. **Generate a sample sheet** with `viralunity create-samplesheet`
2. **Download databases** (for meta pipeline) with `viralunity get-databases`
3. **Run the pipeline** with `viralunity consensus` or `viralunity meta`

```{tip}
Always use absolute paths to avoid mistakes when specifying file locations.
```
