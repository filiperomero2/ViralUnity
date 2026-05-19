# Commands Reference

## `viralunity create-samplesheet`

Before running any pipeline you need a CSV sample sheet. The `create-samplesheet` command generates it automatically from a sequencing run directory.

```bash
viralunity create-samplesheet --input /path/to/run/ --output samples.csv
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--input` | *(required)* | Run directory: contains FASTQ files directly (level 0) or one subdirectory per sample (level 1). |
| `--output` | *(required)* | Output CSV file path. |
| `--level` | `1` | `0` = files in `--input`, `1` = files in subdirectories (one per sample). |
| `--separator` | `-` | Character used to split the file/directory name and extract the sample ID (`-`, `_`, or `.`). |
| `--pattern` | `R1` | Pattern identifying the first read file at level 0 (`R1` for Illumina, `barcode` for Nanopore). |

### Output format

The generated CSV has no header row:

```
# Illumina (3 columns)
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz

# Nanopore (2 columns)
barcode01,/path/to/barcode01.fastq
```

### Examples

**Illumina** — FASTQs organised in one subdirectory per sample (level 1, default):

```bash
viralunity create-samplesheet \
    --input /path/to/illumina_run/ \
    --output /path/to/samples.csv
```

**Illumina** — all FASTQs directly in the run directory (level 0):

```bash
viralunity create-samplesheet \
    --input /path/to/illumina_run/ \
    --output /path/to/samples.csv \
    --level 0 \
    --separator _ \
    --pattern R1
```

**Nanopore** — one FASTQ per barcode subdirectory:

```bash
viralunity create-samplesheet \
    --input /path/to/nanopore_run/ \
    --output /path/to/samples_nano.csv \
    --separator _ \
    --pattern barcode
```

---

## `viralunity consensus`

The consensus pipeline takes raw reads to processed consensus genome sequences with a single command. Select the data type as a subcommand (`illumina` or `nanopore`).

### Options — shared (both data types)

| Option | Default | Description |
|--------|---------|-------------|
| `--sample-sheet` | *(required)* | CSV file with sample IDs and file paths. |
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

### Options — Illumina only

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

### Options — Nanopore only

| Option | Default | Description |
|--------|---------|-------------|
| `--chunk-size` | `10000` | Chunk size for clair3 processing. |
| `--clair3-model` | `r1041_e82_400bps_sup_v500` | Clair3 model for variant calling. |
| `--variant-quality` | `20` | Minimum variant quality (clair3). |
| `--variant-depth` | `10` | Minimum alt allele depth (clair3). |
| `--minimum-map-quality` | `30` | Minimum mapping quality (clair3). |

### Examples

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

---

## `viralunity meta`

The metagenomics pipeline takes raw reads to taxonomic classifications and visualizations. You can use **Kraken2 only**, **Diamond only**, or **both**. Each tool can be run on **reads** and, when assembly is enabled, on **contigs**.

### Pipeline overview (Illumina)

1. **Quality control** — [fastp](https://github.com/OpenGene/fastp): trimming, quality filtering, adapter handling, JSON/HTML reports.
2. **Optional dehosting** — minimap2 (if `--host-reference` is set) or [Deacon](https://github.com/bede/deacon) (if `--deacon-index` is set).
3. **Read classification** — [Kraken2](https://github.com/DerrickWood/kraken2) and/or [DIAMOND](https://github.com/bbuchfink/diamond) blastx on reads (optional each).
4. **Optional de novo assembly** — [MEGAHIT](https://github.com/voutcn/megahit) on host-filtered pairs.
5. **Contig classification** — Kraken2 and/or Diamond on contigs (when assembly is run).
6. **Summaries and filters** — Per-sample taxa tables, RPM normalization, optional max-RPM bleed filter and negative-control background filter. [Krona](https://github.com/marbl/Krona) plots and [MultiQC](https://multiqc.info) reports.
7. **Optional dynamic reference assembly** — Automatic reference sequence selection from taxonomic hits and subsequent consensus assembly.

### Pipeline overview (Nanopore)

1. **Optional dehosting** — minimap2 or Deacon.
2. **Read classification** — Kraken2 and/or Diamond on reads.
3. **Optional de novo assembly** — MEGAHIT on host-filtered reads.
4. **Optional polishing** — [Racon](https://github.com/lbcb-sci/racon) and/or [Medaka](https://github.com/nanoporetech/medaka) on MEGAHIT contigs.
5. **Contig classification** — Kraken2 and/or Diamond on contigs.
6. **Summaries and filters** — Same as Illumina.
7. **Optional dynamic reference assembly** — Same as Illumina.

### Options — shared (both data types)

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
| `--run-reference-assembly`/`--no-run-reference-assembly` | off | Enable reference assembly from filtered taxonomic hits. |
| `--method` | `kraken2` | Method used for reference assembly (`kraken2`, `diamond`, `both`). Required when `--run-reference-assembly` is set. |
| `--source` | `reads` | Source of taxonomy data for reference assembly (`reads`, `contigs`, `both`). Required when `--run-reference-assembly` is set. |
| `--reads-count` | `100` | Minimum reads assigned to a viral family to trigger reference assembly. |
| `--contigs-count` | `1` | Minimum contigs assigned to a viral family to trigger reference assembly. |
| `--families` | `Coronaviridae,...` | Comma-separated list of targeted viral families for reference assembly. |
| `--reference-selection-strategy` | `taxid` | Strategy to associate a reference genome (`taxid` or `similarity`). |
| `--blast-qcov` | `80` | Minimum query coverage for similarity strategy. |
| `--blast-pident` | `80` | Minimum percent identity for similarity strategy. |
| `--viral-genomes` | `NA` | FASTA of viral genomes for reference assembly database. |
| `--viral-taxids` | `NA` | TSV of genome to TaxID mapping for reference assembly. |
| `--threads` | `1` | Threads per task. |
| `--threads-total` | `1` | Total threads for the workflow. |
| `--create-config-only` | off | Only generate the config; do not run the workflow. |

### Options — Illumina only (fastp QC)

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

### Options — Nanopore only (polishing)

| Option | Default | Description |
|--------|---------|-------------|
| `--run-polish-racon/--no-polish-racon` | off | Racon polishing on MEGAHIT assembly. |
| `--run-polish-medaka/--no-polish-medaka` | off | Medaka polishing on assembly. |
| `--medaka-model` | — | Medaka model (e.g. `r941_min_high_g360`). |

### Examples

**Illumina — Kraken2 on reads (default):**

```bash
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

**Illumina — with reference assembly (Kraken2 hits drive consensus assembly):**

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --run-reference-assembly \
    --method kraken2 \
    --source reads \
    --families Coronaviridae,Orthomyxoviridae \
    --reads-count 100 \
    --viral-genomes databases/virus_genomes/viral.genomes.fasta \
    --viral-taxids databases/virus_genomes/genome2taxid.tsv \
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

---

## `viralunity get-databases`

The `get-databases` command downloads and sets up the reference databases required by the meta pipeline. Each subcommand takes an optional `--path` argument (default: `databases` in the current directory), inside which a named subdirectory is created.

```bash
viralunity get-databases --help
viralunity get-databases kraken2 --help
```

### `kraken2` — Kraken2 pre-built index

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/kraken2/`. |
| `--url` | k2_viral_20240112 | URL of the pre-built index archive. Check [available indexes](https://benlangmead.github.io/aws-indexes/k2). |

```bash
viralunity get-databases kraken2 --path /data/dbs
# database at /data/dbs/kraken2/
# use: --kraken2-database /data/dbs/kraken2
```

### `krona` — Krona taxonomy

Requires the viralunity conda environment to be active (`CONDA_PREFIX` must be set).

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/krona/taxonomy/`. |

```bash
conda activate viralunity
viralunity get-databases krona --path /data/dbs
# taxonomy at /data/dbs/krona/taxonomy/
# use: --krona-database /data/dbs/krona/taxonomy
```

### `taxdump` — NCBI taxdump

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/taxdump/`. |
| `--url` | NCBI taxdump | URL of the taxdump archive. |

```bash
viralunity get-databases taxdump --path /data/dbs
# nodes.dmp, names.dmp, etc. at /data/dbs/taxdump/
# use: --taxdump /data/dbs/taxdump
```

### `diamond` — Diamond protein database

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/diamond/`. |
| `--taxon` | Viruses | NCBI taxon name to download (e.g. `Viruses`, `coronaviridae`). |
| `--refseq/--no-refseq` | `on` | Limit to RefSeq genomes only. |
| `--threads` | `1` | Threads for `diamond makedb`. |
| `--skip-makedb` | off | Download and reformat files only; skip `diamond makedb`. |

Requires the NCBI Datasets CLI (`conda install -c conda-forge ncbi-datasets-cli`).

```bash
viralunity get-databases diamond --path /data/dbs --threads 4
# /data/dbs/diamond/viral.dmnd
# /data/dbs/diamond/protein2taxid.tsv
# use: --diamond-database /data/dbs/diamond/viral.dmnd
#      --taxids /data/dbs/diamond/protein2taxid.tsv
```

NCBI Datasets occasionally bundles a small number of nucleotide CDS records
(with a genome accession such as `NC_xxxxxx.1` as the header) inside the
`protein.faa` it returns for `--include protein`. The `diamond` subcommand
automatically detects and drops these records before running
`diamond makedb` and prints a warning summarizing how many were skipped.

### `clean-protein-fasta` — Strip nucleotide records from a protein FASTA

| Option | Default | Description |
|--------|---------|-------------|
| `--input` | *(required)* | Input protein FASTA file (e.g. `databases/diamond/viral.protein.faa`). |
| `--output` | `<input>.cleaned.faa` | Output FASTA path. If equal to `--input`, the file is replaced atomically. |
| `--no-backup` | off | When replacing the input in place, do not keep a `.with_dna.bak` backup. |
| `--min-dna-length` | `20` | Minimum sequence length to classify a record as DNA. |

Use this command to recover an already-downloaded `viral.protein.faa` that
fails `diamond makedb` with:

```
Error: The sequences are expected to be proteins but only contain DNA letters.
```

```bash
viralunity get-databases clean-protein-fasta \
    --input databases/diamond/viral.protein.faa \
    --output databases/diamond/viral.protein.faa

diamond makedb \
    --in databases/diamond/viral.protein.faa \
    --db databases/diamond/viral.dmnd \
    --threads 4
```

### `virus-genome` — Viral genomes database

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/virus_genomes/`. |
| `--taxon` | Viruses | NCBI taxon name to download (e.g. `Viruses`, `coronaviridae`). |
| `--refseq`/`--no-refseq` | `on` | Limit to RefSeq genomes only. |
| `--skip-makeblastdb` | off | Download and reformat files only; skip `makeblastdb` index creation. |

Requires the NCBI Datasets CLI and BLAST+ (`makeblastdb`).

After downloading and reformatting genome sequences, the command automatically runs `makeblastdb` to build a nucleotide BLAST index alongside the FASTA. This index is required when using the `similarity` reference selection strategy (`--reference-selection-strategy similarity`).

```bash
viralunity get-databases virus-genome --taxon Viruses
# FASTA at databases/virus_genomes/viral.genomes.fasta
# BLAST DB index files at databases/virus_genomes/viral.genomes.fasta.{nhr,nin,nsq,...}
# taxids at databases/virus_genomes/genome2taxid.tsv
# use: --viral-genomes databases/virus_genomes/viral.genomes.fasta
#      --viral-taxids databases/virus_genomes/genome2taxid.tsv
```

### `host-genome` — Host genome download

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/host_genomes/`. |
| `--accession` | *(required)* | NCBI genome accession ID (e.g. `GCA_000001405.29`). |

```bash
viralunity get-databases host-genome --accession GCA_000001405.29
# fasta at databases/host_genomes/GCA_000001405.29.fasta
```

### `deacon-index` — Pre-built Deacon index

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/deacon_indexes/`. |
| `--index-name` | `panhuman-1` | Pre-built index name (`panhuman-1` or `panmouse-1`). |

```bash
viralunity get-databases deacon-index --index-name panhuman-1
# index at databases/deacon_indexes/panhuman-1.idx
# use: --deacon-index databases/deacon_indexes/panhuman-1.idx
```

### Download all databases at once

```bash
viralunity get-databases all --threads 4
```

---

## `viralunity build-deacon-index`

Build a Deacon minimizer index from a host genome FASTA.

| Option | Default | Description |
|--------|---------|-------------|
| `--path` | `databases` | Parent directory; creates `{path}/deacon_indexes/`. |
| `--input` | *(required)* | Input FASTA file to index. |
| `--threads` | `8` | Number of threads to use. |

```bash
viralunity build-deacon-index \
    --input databases/host_genomes/GCA_000001405.29.fasta \
    --threads 4
# output at databases/deacon_indexes/GCA_000001405.29.idx
```
