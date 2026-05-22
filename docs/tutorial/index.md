# Tutorial

This tutorial walks you from a fresh ViralUnity install to a finished consensus genome or a finished metagenomic profile on real example data. It assumes no prior experience with the tool. By the end, you will know how to point either pipeline at your own viral high-throughput sequencing data and interpret the output.

ViralUnity orchestrates a Snakemake workflow that turns raw viral sequencing reads into either:

- **a polished consensus genome per sample**, when you already know what you are sequencing and have a reference (`viralunity consensus`); or
- **a taxonomic profile** of what is in each sample — plus, optionally, de novo contigs and reference-guided consensus genomes for any virus that is detected — when the contents are unknown or mixed (`viralunity meta`).

## Pick a pipeline

| If you have…                                                                  | …use                                              |
|-------------------------------------------------------------------------------|---------------------------------------------------|
| Reads from a known virus and a reference genome; want a per-sample consensus  | `viralunity consensus`                            |
| Reads from a clinical/environmental sample; want to know what viruses are in it | `viralunity meta` (defaults)                      |
| Same as above, but also want consensus genomes for any virus that turns up    | `viralunity meta --run-reference-assembly`        |
| Multi-segment virus (influenza, bunyaviruses, …)                              | `viralunity consensus --segmented-reference …`    |

## How to read this tutorial

1. [Setup](setup.md) — install ViralUnity, generate sample sheets, and download the reference databases needed by the metagenomic pipeline. **Do this once.**
2. [Consensus pipeline](consensus.md) — walk through reference-guided consensus assembly on Illumina and Nanopore SARS-CoV-2 data, plus how to handle segmented viruses.
3. [Metagenomics pipeline](metagenomics.md) — build up the metagenomic workflow incrementally, from a minimal Kraken2 classification to the full dehosting + assembly + DIAMOND + reference-assembly pipeline. Includes a deep dive into the contamination filters and the dynamic reference-assembly step.

Each page is self-contained but assumes you have completed [Setup](setup.md) first.

```{toctree}
:maxdepth: 1
:hidden:

setup
consensus
metagenomics
```

## About the example data

The worked examples use a small SARS-CoV-2 dataset that ships in this checkout under `my_test_data/`:

- `my_test_data/illumina_data/` — two paired-end Illumina amplicon samples: `4117_S80_L001_R{1,2}_001.fastq.gz` and `61382_S1_L001_R{1,2}_001.fastq.gz`.
- `my_test_data/nanopore_data/` — two Nanopore samples: `barcode05.fastq` and `barcode09.fastq`.

```{note}
`my_test_data/` is gitignored and only present in checkouts where you have placed it yourself. If you do not have it, substitute your own SARS-CoV-2 FASTQs in the commands below — the paths and sample IDs are the only things that change. A suitable reference genome is NCBI [`MN908947.3`](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3), available as `nCoV-2019.reference.fasta` in any of the public [ARTIC](https://github.com/artic-network/primer-schemes) primer-scheme repositories.
```

For matching example outputs, the repository also includes `my_results/test_consensus_*/` and `my_results/test_meta_*/` — these are the outputs produced by running each tutorial example. Use them to sanity-check what your own run should look like.

## Further reading

The tutorial covers the common path through ViralUnity. Once you are productive, the reference pages have the exhaustive detail:

- [Commands reference](../commands.md) — every CLI option, with defaults and units.
- [Output layout](../output.md) — the full directory tree both pipelines produce.
- [Notes](../notes.md) — segmented viruses, the two reference-selection strategies for dynamic reference assembly, and other advanced topics.
