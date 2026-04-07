# ViralUnity

ViralUnity is a simple tool to perform analysis of viral high-throughput sequencing data. It comprises a collection of Python scripts and Snakemake workflows, including several software dependencies, to perform data quality control, taxonomic assignments, and reference genome assembly.

ViralUnity runs on *nix systems and is able to process entire sequencing runs in minimal time on a regular computer.

## Installation

ViralUnity is a Python package that launches Snakemake workflows. The tool dependencies are documented in the conda environment file `environment.yml` and the workflows dependencies are documented in the `viralunity/scripts/envs` directory.

To enable ViralUnity, clone the repo and create the environment:

```bash
git clone https://github.com/InstitutoTodosPelaSaude/ViralUnity.git
cd ViralUnity/
conda env create -f environment.yml
conda activate viralunity
```

Or with micromamba (recommended on macOS with Apple Silicon):

```bash
micromamba env create -f environment.yml --platform osx-64
```

> [!WARNING]
> On macOS with Apple Silicon (M1 or later), the current `viralunity/scripts/envs/clair3.yaml` may not work due to clair3 dependencies.

## Usage

ViralUnity provides five main subcommands:

```
viralunity [--version]
├── build-deacon-index   Build a Deacon minimizer index from a FASTA file
├── create-samplesheet   Generate a sample-sheet CSV from a run directory
├── get-databases        Download and set up reference databases
│   ├── kraken2          Download a Kraken2 pre-built index
│   ├── krona            Set up the Krona taxonomy database
│   ├── taxdump          Download the NCBI taxdump
│   ├── diamond          Download RefSeq viral proteins and build a Diamond DB
│   ├── host-genome      Download a host genome using NCBI Datasets
│   └── deacon-index     Download a pre-built Deacon minimizer index
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

---

### Generate a sample sheet (`viralunity create-samplesheet`)

Before running any pipeline you need a CSV sample sheet. The `create-samplesheet` command generates it automatically from a sequencing run directory.

```bash
viralunity create-samplesheet --input /path/to/run/ --output samples.csv
```

#### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--input` | *(required)* | Run directory: contains FASTQ files directly (level 0) or one subdirectory per sample (level 1). |
| `--output` | *(required)* | Output CSV file path. |
| `--level` | `1` | `0` = files in `--input`, `1` = files in subdirectories (one per sample). |
| `--separator` | `-` | Character used to split the file/directory name and extract the sample ID (`-`, `_`, or `.`). |
| `--pattern` | `R1` | Pattern identifying the first read file at level 0 (`R1` for Illumina, `barcode` for Nanopore). |

#### Output format

The generated CSV has no header row:

```
# Illumina (3 columns)
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz

# Nanopore (2 columns)
barcode01,/path/to/barcode01.fastq
```

#### Examples

**Illumina run** — FASTQs organised in one subdirectory per sample (level 1, default):

```bash
viralunity create-samplesheet \
    --input /path/to/illumina_run/ \
    --output /path/to/samples.csv
```

**Illumina run** — all FASTQs directly in the run directory (level 0):

```bash
viralunity create-samplesheet \
    --input /path/to/illumina_run/ \
    --output /path/to/samples.csv \
    --level 0 \
    --separator _ \
    --pattern R1
```

**Nanopore run** — one FASTQ per barcode subdirectory:

```bash
viralunity create-samplesheet \
    --input /path/to/nanopore_run/ \
    --output /path/to/samples_nano.csv \
    --separator _ \
    --pattern barcode
```

---

### Reference viral genome assembly (`viralunity consensus`)

The consensus pipeline takes raw reads to processed consensus genome sequences with a single command. Select the data type as a subcommand (`illumina` or `nanopore`).

#### Options — shared (both data types)

| Option | Default | Description |
|--------|---------|-------------|
| `--sample-sheet` | *(required)* | CSV file with sample IDs and file paths.  |
| `--config-file` | *(required)* | Path for the YAML config file to be created. |
| `--output` | *(required)* | Base output directory. |
| `--reference` | — | Reference genome FASTA (mutually exclusive with `--segmented-reference`). |
| `--segmented-reference` | — | Per-segment reference: `SEGMENT=PATH` (repeatable). |
| `--primer-scheme` | — | Primer scheme BED file (amplicon sequencing only). |
| `--minimum-coverage` | `20` | Minimum depth for consensus base inclusion. |
| `--minimum-read-length` | `50` | Minimum read length threshold. |
| `--af-threshold` | `0.51` | Min allele frequency to call variant into consensus. |
| `--run-name` | `undefined` | Name for the sequencing run. |
| `--threads` | `1` | Threads per individual task. |
| `--threads-total` | `1` | Total threads for the workflow. |
| `--create-config-only` | off | Only generate the config file; do not run the workflow. |

#### Options — Illumina only

| Option | Default | Description |
|--------|---------|-------------|
| `--adapters` | — | Adapter sequences FASTA (fastp QC). |
| `--trim-head` | `0` | Bases to trim from 5′ end. |
| `--trim-tail` | `0` | Bases to trim from 3′ end. |
| `--cut-front-mean-quality` | `10` | fastp cut_front quality threshold. |
| `--cut-tail-mean-quality` | `10` | fastp cut_tail quality threshold. |
| `--cut-right-window-size` | `4` | fastp cut_right window size. |
| `--cut-right-mean-quality` | `15` | fastp cut_right quality threshold. |
| `--af-isnv-threshold` | `0.0` | Min allele frequency for iSNV analysis. |
| `--run-isnv` | off | Run intra-host SNV analysis (LoFreq). |

#### Options — Nanopore only

| Option | Default | Description |
|--------|---------|-------------|
| `--chunk-size` | `10000` | Chunk size for clair3 processing. |
| `--clair3-model` | `r1041_e82_400bps_sup_v500` | Clair3 model for variant calling. |
| `--variant-quality` | `20` | Minimum variant quality (clair3). |
| `--variant-depth` | `10` | Minimum alt allele depth (clair3). |
| `--minimum-map-quality` | `30` | Minimum mapping quality (clair3). |

#### Run examples

```bash
viralunity create-samplesheet \
    --input /path/to/fastq_dir/ \
    --output /path/to/example.csv
```

**Illumina — single reference:**

```bash
viralunity consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --reference /path/to/references/viral_genome.fasta \
    --primer-scheme /path/to/primers.bed \
    --adapters /path/to/adapters.fa \
    --threads 2 \
    --threads-total 4
```

**Illumina — segmented genome** (replace `--reference` with `--segmented-reference`):

```bash
viralunity consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example_segmented.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --segmented-reference L=/path/to/L_segment.fasta \
    --segmented-reference S=/path/to/S_segment.fasta \
    --threads 2 \
    --threads-total 4
```

**Nanopore:**

```bash
viralunity consensus nanopore \
    --sample-sheet /path/to/samplesheet_nano.csv \
    --config-file /path/to/example_nano.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --reference /path/to/reference.fasta \
    --threads 4 \
    --threads-total 4
```

**Config only** (generates config without running the workflow):

```bash
viralunity consensus illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example.yml \
    --output /path/to/example_output \
    --reference /path/to/reference.fasta \
    --create-config-only
```

The output directory contains QC reports, mapping files, coverage reports, and consensus sequences.

> [!TIP]
> Always use absolute paths to avoid mistakes.

---

### Taxonomic assignment (`viralunity meta`)

The metagenomics pipeline takes raw reads to taxonomic classifications and visualizations. You can use **Kraken2 only**, **Diamond only**, or **both**. Each tool can be run on **reads** and, when assembly is enabled, on **contigs**.

#### Sample sheet format

The `--sample-sheet` must point to a CSV file (no header row):

- **Illumina:** three columns — sample name, R1 path, R2 path.
- **Nanopore:** two columns — sample name, path to single FASTQ/FASTA.

```
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

```
barcode01,/path/to/barcode01.fastq
barcode02,/path/to/barcode02.fastq.gz
```

#### Pipeline overview (Illumina)

1. **Quality control** — [fastp](https://github.com/OpenGene/fastp): trimming, quality filtering, adapter handling, JSON/HTML reports.
2. **Optional dehosting** — If `--host-reference` is set, [minimap2](https://github.com/lh3/minimap2) aligns reads; unmapped reads are kept. Use `--deacon-index` for [Deacon](https://github.com/bede/deacon)-based host depletion (faster, no indexing).
3. **Read classification** — [Kraken2](https://github.com/DerrickWood/kraken2) and/or [DIAMOND](https://github.com/bbuchfink/diamond) blastx on reads (optional each).
4. **Optional de novo assembly** — [MEGAHIT](https://github.com/voutcn/megahit) on host-filtered pairs.
5. **Contig classification** — Kraken2 and/or Diamond on contigs (when assembly is run).
6. **Summaries and filters** — Per-sample taxa tables, RPM normalization, optional max-RPM bleed filter and negative-control background filter. [Krona](https://github.com/marbl/Krona) plots and [MultiQC](https://multiqc.info) reports.

#### Pipeline overview (Nanopore)

1. **Optional dehosting** — minimap2 or Deacon.
2. **Read classification** — Kraken2 and/or Diamond on reads.
3. **Optional de novo assembly** — MEGAHIT on host-filtered reads.
4. **Optional polishing** — [Racon](https://github.com/lbcb-sci/racon) and/or [Medaka](https://github.com/nanoporetech/medaka) on MEGAHIT contigs.
5. **Contig classification** — Kraken2 and/or Diamond on contigs.
6. **Summaries and filters** — Same as Illumina.

#### Options — shared (both data types)

| Option | Default | Description |
|--------|---------|-------------|
| `--sample-sheet` | *(required)* | CSV with sample IDs and file paths. |
| `--config-file` | *(required)* | Path to YAML config file to be generated. |
| `--output` | *(required)* | Base output directory. |
| `--run-name` | `undefined` | Run name; output goes under `{output}/{run_name}/`. |
| `--kraken2-database` | `NA` | Path to Kraken2 database. |
| `--krona-database` | `NA` | Path to Krona taxonomy (required for any classification). |
| `--taxdump` | `NA` | NCBI taxdump dir (`nodes.dmp`, `names.dmp`). |
| `--remove-human-reads` | off | Remove human reads from summaries/Krona. |
| `--remove-unclassified-reads` | off | Remove unclassified reads from summaries/Krona. |
| `--host-reference` | `NA` | Host genome FASTA for minimap2 dehosting. |
| `--deacon-index` | `NA` | Deacon minimizer index for host depletion. |
| `--run-denovo-assembly` | off | Run MEGAHIT de novo assembly. |
| `--run-kraken2-reads/--no-kraken2-reads` | on | Kraken2 on reads. |
| `--run-kraken2-contigs/--no-kraken2-contigs` | on | Kraken2 on contigs. |
| `--run-diamond-reads/--no-diamond-reads` | off | Diamond on reads. |
| `--run-diamond-contigs/--no-diamond-contigs` | off | Diamond on contigs. |
| `--diamond-database` | `NA` | Protein FASTA for Diamond. |
| `--taxids` | `NA` | NCBI taxid mapping for Diamond taxonomy. |
| `--diamond-sensitivity` | `sensitive` | `sensitive` / `mid-sensitive` / `more-sensitive` / `ultra-sensitive`. |
| `--evalue` | `0.001` | Diamond E-value threshold. |
| `--bleed-fraction` | `0.005` | Max-RPM bleed filter fraction. |
| `--negative-controls` | (empty) | Comma-separated sample IDs used as negative controls. |
| `--negative-p-threshold` | `0.01` | p-value threshold for negative-control filter. |
| `--minimum-hit-group` | `4` | Kraken2 minimum-hit-group parameter. |
| `--threads` | `1` | Threads per task. |
| `--threads-total` | `1` | Total threads for the workflow. |
| `--create-config-only` | off | Only generate the config; do not run the workflow. |

#### Options — Illumina only (fastp QC)

| Option | Default | Description |
|--------|---------|-------------|
| `--adapters` | `NA` | Adapter FASTA; `NA` = auto-detect. |
| `--minimum-read-length` | `50` | Minimum read length after trimming. |
| `--trim-head` | `0` | Bases to trim from 5′. |
| `--trim-tail` | `0` | Bases to trim from 3′. |
| `--cut-front-mean-quality` | `20` | fastp cut_front quality threshold. |
| `--cut-tail-mean-quality` | `20` | fastp cut_tail quality threshold. |
| `--cut-right-window-size` | `4` | fastp cut_right window size. |
| `--cut-right-mean-quality` | `20` | fastp cut_right quality threshold. |

#### Options — Nanopore only (polishing)

| Option | Default | Description |
|--------|---------|-------------|
| `--run-polish-racon/--no-polish-racon` | off | Racon polishing on MEGAHIT assembly. |
| `--run-polish-medaka/--no-polish-medaka` | off | Medaka polishing on assembly. |
| `--medaka-model` | — | Medaka model (e.g. `r941_min_high_g360`). |

### Download databases (`viralunity get-databases`)

Some steps into `meta` command requires reference databases, others require host genome or indexes to filter host reads.

The `get-databases` command downloads and sets up the reference databases required by workflows. Each subcommand takes an optional `--path` argument (default: `databases` in the current directory) inside which a named subdirectory is created (e.g. `--path /data/dbs` creates `/data/dbs/kraken2/`, `/data/dbs/krona/`, etc.).

```bash
viralunity get-databases --help
viralunity get-databases kraken2 --help
```

#### `kraken2` — Kraken2 pre-built index

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/kraken2/`. |
| `--url` | k2_viral_20240112 | URL of the pre-built index archive. Check [available indexes](https://benlangmead.github.io/aws-indexes/k2). |

```bash
viralunity get-databases kraken2 --path /data/dbs
# database at /data/dbs/kraken2/
# use: --kraken2-database /data/dbs/kraken2
```

#### `krona` — Krona taxonomy

Requires the viralunity conda environment to be active (`CONDA_PREFIX` must be set).

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/krona/taxonomy/`. |

```bash
mamba activate viralunity
viralunity get-databases krona --path /data/dbs
# taxonomy at /data/dbs/krona/taxonomy/
# use: --krona-database /data/dbs/krona/taxonomy
```

The command removes the stale taxonomy bundled with the conda package, symlinks the new directory into `$CONDA_PREFIX/opt/krona/taxonomy`, and runs `ktUpdateTaxonomy.sh`.

#### `taxdump` — NCBI taxdump

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/taxdump/`. |
| `--url` | NCBI taxdump | URL of the taxdump archive. |

```bash
viralunity get-databases taxdump --path /data/dbs
# nodes.dmp, names.dmp, etc. at /data/dbs/taxdump/
# use: --taxdump /data/dbs/taxdump
```

#### `diamond` — Diamond protein database

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/diamond/`. |
| `--taxon` | Viruses | NCBI taxon name to download (e.g. 'Viruses', 'coronaviridae'). |
| `--refseq/--no-refseq` | `on` | Limit to RefSeq genomes only. |
| `--threads` | `1` | Threads for `diamond makedb`. |
| `--skip-makedb` | off | Download and reformat files only; skip `diamond makedb`. |

Requires the NCBI Datasets CLI to be installed (`conda install -c conda-forge ncbi-datasets-cli`). The default database is built with RefSeq genomes only. To include all viral genomes, use the `--no-refseq` flag.

```bash
viralunity get-databases diamond --path /data/dbs --threads 4
# /data/dbs/diamond/viral.dmnd
# /data/dbs/diamond/protein2taxid.tsv
# use: --diamond-database /data/dbs/diamond/viral.dmnd
#      --taxids /data/dbs/diamond/protein2taxid.tsv
```

#### Download all databases at once

Instead of running each get-databases command individually, you can run the `all` command to download all databases at once.

```bash
viralunity get-databases all --threads 4
```

#### Filter databases

The host filter features allow you to filter out host reads from your samples. This is useful for metagenomic analysis where you want to focus on viral reads.

For user convenience we included commands to download genomes from NCBI as well as deacon pre-built indexes, but you can also build or provide your own host files.

##### `host-genome` — Host Genome

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/host_genomes/`. |
| `--accession` | *(required)* | NCBI genome accession ID to download (e.g. `GCA_000001405.29`). |

Requires the NCBI Datasets CLI to be installed.

```bash
viralunity get-databases host-genome --accession GCA_000001405.29
# fasta at databases/host_genomes/GCA_000001405.29.fasta
# info at databases/host_genomes/GCA_000001405.29.info.txt
```

##### `deacon-index` — Pre-built Deacon Index

For user convenience, we include the command `deacon index fetch` inside viralunity, but feel free to use deacon commands themselves if you prefer.

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/deacon_indexes/`. |
| `--index-name` | `panhuman-1` | Pre-built index name to download (`panhuman-1` or `panmouse-1`). |


```bash
viralunity get-databases deacon-index --index-name panhuman-1
# index at databases/deacon_indexes/panhuman-1.idx
# use: --deacon-index databases/deacon_indexes/panhuman-1.idx
```

---

### Build a Deacon Index (`viralunity build-deacon-index`)

If you downloaded a host genome FASTA, you can build a Deacon index from it. 

For user convenience, we include the command `deacon index build` inside viralunity, but feel free to use deacon commands themselves if you prefer.

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/deacon_indexes/`. |
| `--input` | *(required)* | Input FASTA file to index. |
| `--threads` | `8` | Number of threads to use. |

```bash
viralunity build-deacon-index --input databases/host_genomes/GCA_000001405.29.fasta --threads 4
# output at databases/deacon_indexes/GCA_000001405.29.idx
```

---


#### Run examples

**Illumina — Kraken2 on reads (default):**

```bash
conda activate viralunity
viralunity meta illumina \
    --sample-sheet /path/to/example.csv \
    --config-file /path/to/example_meta.yml \
    --run-name example_run \
    --output /path/to/example_output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --threads 2 \
    --threads-total 4
```

**Illumina — QC only** (no classification):

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --no-kraken2-reads
```

**Illumina — with adapters and dehosting:**

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --adapters /path/to/adapters.fa \
    --host-reference /path/to/host_genome.fa \
    --threads 4 --threads-total 8
```

**Illumina — Diamond on reads only** (no Kraken2):

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --no-kraken2-reads \
    --run-diamond-reads \
    --diamond-database /path/to/proteins.faa \
    --taxids /path/to/protein2taxid.tsv \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump
```

**Illumina — full pipeline (assembly + Kraken2 + Diamond, reads and contigs):**

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --run-diamond-reads --run-diamond-contigs \
    --diamond-database /path/to/proteins.faa \
    --taxids /path/to/protein2taxid.tsv \
    --run-denovo-assembly \
    --host-reference /path/to/host.fa \
    --threads 4 --threads-total 8
```

**Nanopore — full pipeline (dehosting, assembly, Medaka polishing, Kraken2 + Diamond):**

```bash
viralunity meta nanopore \
    --sample-sheet /path/to/samplesheet_nano.csv \
    --config-file /path/to/config_nano.yml \
    --run-name example_run \
    --output /path/to/output_nano \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --diamond-database /path/to/proteins.faa \
    --taxids /path/to/protein2taxid.tsv \
    --host-reference /path/to/host.fa \
    --run-denovo-assembly \
    --run-kraken2-reads --run-kraken2-contigs \
    --run-diamond-reads --run-diamond-contigs \
    --run-polish-medaka --medaka-model r941_min_high_g360 \
    --threads 4 --threads-total 4
```

Output includes `qc/` (fastp reports + FastQC + MultiQC — Illumina only), `host_filtered/` (if dehosting is used), `metagenomics/taxonomic_assignments/` (Kraken2 and/or Diamond results, Krona plots, taxa summaries), and optionally `denovo_assembly/` and `medaka_work/` for Nanopore polishing.

---

## Notes

### Segmented viruses

ViralUnity natively supports the assembly of segmented viral genomes. Instead of running the pipeline multiple times, pass multiple segment references with `--segmented-reference`:

```bash
viralunity consensus illumina \
    --sample-sheet samples.csv \
    --config-file config_segmented.yml \
    --output /path/to/output \
    --segmented-reference S=/path/to/s_segment.fasta \
    --segmented-reference L=/path/to/l_segment.fasta
```

When this argument is detected, the pipeline dispatches a specialized modular workflow that processes each segment independently in parallel. Results are organized under `samples/{sample_name}/{segment_name}/`.

### Nanopore data

The metagenomics pipeline supports Nanopore data via `viralunity meta nanopore`. There is no fastp-based QC step; optional dehosting, MEGAHIT assembly, and Racon/Medaka polishing are available. Sample sheets have two columns (sample ID, path to one FASTQ/FASTA per sample).

---

## Tests

```bash
pytest -v test
```

---

## Citation

A scientific publication fully describing this pipeline is being prepared. Meanwhile, feel free to cite this GitHub repository. Primary references for the dependencies used should also be cited:

<a href="https://github.com/OpenGene/fastp">fastp</a>: Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018; 34(17):i884–i890.

<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">FastQC</a>

<a href="https://github.com/ewels/MultiQC">MultiQC</a>: Ewels P, Magnusson M, Lundin S, et al. MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016; 32(19):3047–8.

<a href="https://github.com/lh3/minimap2">Minimap2</a>: Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018; 34:3094-3100.

<a href="https://github.com/samtools/samtools">Samtools</a>: Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009; 25(16):2078–9.

<a href="https://github.com/andersen-lab/ivar">iVar</a>: Grubaugh ND, Gangavarapu K, Quick J, et al. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biol 20, 8 (2019).

<a href="https://github.com/arq5x/bedtools2">BEDtools</a>: Quinlan AR, Hall IM. BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics. 2010; 26(6):841–2.

<a href="https://doi.org/10.1186/s13059-019-1891-0">Kraken2</a>: Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019).

<a href="https://github.com/bbuchfink/diamond">DIAMOND</a>: Buchfink B, Reuter K, Drost H-G. Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods 18, 366–368 (2021).

<a href="https://github.com/voutcn/megahit">MEGAHIT</a>: Li D, Liu C-M, Luo R, Sadakane K, Lam T-W. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics. 2015; 31(10):1674–1676.

<a href="https://doi.org/10.1186/1471-2105-12-385">Krona</a>: Ondov BD, Bergman NH, Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics 12, 385 (2011).

For visualization of FASTA and BAM files, we recommend <a href="https://ormbunkar.se/aliview/">AliView</a> and <a href="https://ics.hutton.ac.uk/tablet/">Tablet</a>, respectively.
