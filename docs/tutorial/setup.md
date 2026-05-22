# Setup

This page covers everything you need before running either pipeline: installing ViralUnity, generating sample sheets for the example data, and downloading the reference databases. Most of it is one-time work; once you finish, you can move straight to either pipeline walkthrough.

## 1. Install ViralUnity

ViralUnity is distributed alongside a conda environment that pins all of the underlying bioinformatics tools (Snakemake, minimap2, fastp, Kraken2, DIAMOND, Krona, clair3, Medaka, …). The per-rule conda environments under `viralunity/scripts/envs/` are picked up automatically by Snakemake at runtime.

```bash
git clone https://github.com/InstitutoTodosPelaSaude/ViralUnity.git
cd ViralUnity
conda env create -n viralunity -f environment.yml
conda activate viralunity
pip install -e .
```

Verify:

```bash
viralunity --version
viralunity --help
```

See [Installation](../installation.md) for the macOS Apple Silicon caveat (the `clair3` environment can be sensitive there) and for development installs (`pip install -e ".[dev]"`).

## 2. Generate sample sheets

A *sample sheet* is a no-header CSV that tells ViralUnity which FASTQ files belong to which sample:

- **Illumina** — 3 columns: `sample_id,R1_path,R2_path`
- **Nanopore** — 2 columns: `sample_id,fastq_path`

You can write the CSV by hand, but `viralunity create-samplesheet` builds it for you by scanning a run directory.

The bundled test data has all FASTQs directly inside the run directory (rather than one subdirectory per sample), so we pass `--level 0` and a small naming convention:

```bash
# Illumina: FASTQs in my_test_data/illumina_data/, sample ID = portion before first '_'
viralunity create-samplesheet \
    --input  my_test_data/illumina_data/ \
    --output samples_illumina.csv \
    --level 0 \
    --separator _ \
    --pattern R1
```

This writes:

```text
4117,my_test_data/illumina_data/4117_S80_L001_R1_001.fastq.gz,my_test_data/illumina_data/4117_S80_L001_R2_001.fastq.gz
61382,my_test_data/illumina_data/61382_S1_L001_R1_001.fastq.gz,my_test_data/illumina_data/61382_S1_L001_R2_001.fastq.gz
```

```bash
# Nanopore: FASTQs in my_test_data/nanopore_data/, sample ID = portion before first '.'
viralunity create-samplesheet \
    --input  my_test_data/nanopore_data/ \
    --output samples_nanopore.csv \
    --level 0 \
    --separator . \
    --pattern barcode
```

Which writes:

```text
barcode05,my_test_data/nanopore_data/barcode05.fastq
barcode09,my_test_data/nanopore_data/barcode09.fastq
```

```{tip}
The reference sample sheets that ship with this checkout are at `my_results/samples_illumina.csv` and `my_results/samples_nanopore.csv` — use them to compare against what you just generated.
```

For your own data, the same command works with `--level 1` (one subdirectory per sample, the default) if your run is organised that way; see the [Commands reference](../commands.md#viralunity-create-samplesheet) for the full table of options.

## 3. Decide which databases you need

Different pipelines need different reference data:

| Pipeline                                              | Required                                                           | Optional (enable extra features)                                                  |
|-------------------------------------------------------|--------------------------------------------------------------------|------------------------------------------------------------------------------------|
| `consensus`                                           | reference FASTA                                                    | primer-scheme BED (amplicon data), adapter FASTA (custom adapters)                |
| `meta` (minimum)                                      | Kraken2 index, Krona taxonomy                                      | NCBI taxdump (filters), host FASTA or Deacon index (dehosting)                    |
| `meta` + DIAMOND                                      | Kraken2, Krona, taxdump, **DIAMOND DB**, **DIAMOND taxids**         | (same as above)                                                                    |
| `meta` + `--run-reference-assembly`                   | the above, **viral genomes FASTA + `genome2taxid.tsv`**, taxdump   | BLAST index (built automatically) — for the `similarity` selection strategy        |

The consensus pipeline does not need any database downloads — only a reference FASTA you already have. The rest of this section sets you up for the metagenomic walkthrough.

## 4. Download the metagenomic databases

The shortcut command grabs the four most-used databases in one go (Kraken2 + Krona + taxdump + DIAMOND). The conda environment must be active because the Krona step needs `ktUpdateTaxonomy.sh`:

```bash
conda activate viralunity
viralunity get-databases all --path databases/ --threads 4
```

Plan for the download. Approximate sizes on disk after extraction:

| Database       | Path created                              | Approx. size | Used by                                                          |
|----------------|-------------------------------------------|--------------|------------------------------------------------------------------|
| Kraken2 viral  | `databases/kraken2/`                      | ~1.1 GB      | `--kraken2-database`                                             |
| Krona taxonomy | `databases/krona/taxonomy/`               | ~150 MB      | `--krona-database`                                               |
| NCBI taxdump   | `databases/taxdump/`                      | ~600 MB      | `--taxdump` (filters, reference assembly)                        |
| DIAMOND viral  | `databases/diamond/`                      | ~400 MB      | `--diamond-database`, `--taxids`                                 |

If you only plan to use Kraken2 on reads, you can skip `get-databases diamond` and only run the first three. Each subcommand is also available standalone — see [`viralunity get-databases --help`](../commands.md#viralunity-get-databases).

### Two more databases for the full meta walkthrough

The `all` subcommand intentionally leaves out two databases that are large and only needed for specific features:

```bash
# Viral genomes + per-genome taxid map + BLAST index — required for --run-reference-assembly
viralunity get-databases virus-genome --path databases/

# Deacon human host index — recommended for dehosting (faster than minimap2 vs a host FASTA)
viralunity get-databases deacon-index --path databases/ --index-name panhuman-1
```

`virus-genome` ships with `viral.genomes.fasta`, `genome2taxid.tsv`, and the BLAST `.nhr`/`.nin`/`.nsq` index files (used for the `similarity` reference-selection strategy). The Deacon index lands at `databases/deacon_indexes/panhuman-1.idx`.

```{tip}
If you want to dehost against a non-human genome, use `viralunity get-databases host-genome --accession <NCBI accession>` to download a host FASTA, then either pass it directly via `--host-reference` (minimap2 dehosting) or build a Deacon index with `viralunity build-deacon-index --input <fasta>` for a faster path.
```

## 5. Verify the setup

After downloading, you should have:

```text
databases/
├── kraken2/           # Kraken2 index files
├── krona/taxonomy/    # taxonomy.tab and friends
├── taxdump/           # nodes.dmp, names.dmp, …
├── diamond/           # viral.dmnd, protein2taxid.tsv
├── virus_genomes/     # viral.genomes.fasta(.nhr/.nin/.nsq), genome2taxid.tsv
└── deacon_indexes/    # panhuman-1.idx
```

A quick sanity check:

```bash
ls databases/
ls databases/diamond/ databases/virus_genomes/
```

You are now ready to:

- [Run the consensus pipeline](consensus.md), or
- [Run the metagenomics pipeline](metagenomics.md).
