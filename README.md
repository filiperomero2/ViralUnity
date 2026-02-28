# ViralUnity

ViralUnity is a simple tool to perform analysis of viral high-throughput sequencing data. It comprises a collection of python scripts and snakemake workflows, including several software dependencies, to perform data quality control, taxonomic assignments and reference genome assembly.

ViralUnity runs on *nix systems and is able to process entire sequencing runs in minimal time on a regular computer.

## Installation

ViralUnity is a Python package that launches Snakemake workflows. All dependencies are documented in the conda environment file `environment.yml`. We recommend using <a href="https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html">mamba</a> for installation.

To enable ViralUnity, clone the repo and create the environment:

    git clone https://github.com/filiperomero2/ViralUnity.git
    cd ViralUnity
    conda create -n viralunity python=3.11
    conda activate viralunity
    conda env update -n viralunity --file environment.yml
    pip install .

Tip: On macOS (with Apple Silicon M1 or later), you may need to further configure the environment before installing dependencies:

    conda create -n viralunity
    conda activate viralunity
    conda config --env --set subdir osx-64
    conda env update -n viralunity --file environment.yml
    pip install .

Verify that all required dependencies (listed in `environment.yml`) are available before running the pipeline. On certain systems it may be necessary to manually set conda channel priority.

## Usage

ViralUnity provides two main pipelines:

| Command | Description |
|--------|-------------|
| `viralunity meta` | Metagenomics: taxonomic assignment of reads and/or contigs (Kraken2, Diamond), optional dehosting, assembly, and summaries. |
| `viralunity consensus` | Reference-based consensus assembly from mapped reads. |

Use `viralunity -h` and `viralunity meta -h` or `viralunity consensus -h` for full option lists. 

### Reference viral genome assembly

The viralunity consensus sequence pipeline was designed to be as simple as possible, with the objective of making users go from raw reads to processed consensus genome sequences with a single command line. The core script is called viralunity_consensus.py and accepts several arguments:

    --data-type             Sequencing data type (illumina or nanopore)
    --sample-sheet          Complete path for a csv file with samples data paths and metadata
    --config-file           Complete path for input (viralunity config) file to be created.
    --output                Complete path for output directory to be created by viralunity.
    --run-name              Name for the sequencing run (optional).
    --reference             Complete path for the reference genome in fasta format
    --primer-scheme         Complete path for the primer scheme bed file (amplicon sequencing only).
    --minimum-coverage      Minimum sequencing coverage for including base in consensus sequence (Default = 20)
    --adapters              Complete path for adapters sequences in fasta format (Illumina QC)
    --minimum-read-length   Minimum read length threshold (Default = 50) (Illumina QC)
    --trim                  Number of bases to trim from the 5' end of reads (Default = 0) (Illumina QC)
    --create-config-only    Only create config file, not running the workflow (boolean)
    --threads               Number of available threads for individual tasks (Default = 1)
    --threads-total         Number of available threads for the entire workflow (Default = 1)
    -h, --help              Show this help message and exit

The samplesheet argument receives the path to a csv file that contains sample names and fastq file paths. This file can be created automatically with the script viralunity_create_samplesheet.py (see its --help option). 

If the specified files paths are correct, the script will generate config file containing all information required to execute the snakemake consensus workflow, allowing optimal parallelization schemes. Beyond creating the config file, the script will execute the workflow, performing data quality control (trimmomatic) (Illumina only), read mapping (minimap2), and consensus sequence inference (samtools). QC reports are generated with fastQC and multiQC. For data generated under amplicon sequencing approaches, the option --primer-schemes takes as argument primer-scheme bed files, which are used for trimming these technical sequences (ivar). The results from all these analysis are stored in the specified output directory. 

#### Run

To run the pipeline, go to the repository directory,activate the conda environment and create a samplesheet file:

    python viralunity/viralunity_create_samplesheet.py --input /home/Desktop/test_data/ --output /home/Desktop/example.csv

Check the contents of the samplesheet file. If the correct paths for fastq files are specified, launch the analysis:

    conda activate viralunity
    viralunity consensus --data-type illumina --sample-sheet /home/Desktop/example.csv --config-file /home/Desktop/example.yml --run-name example_run --output /home/Desktop/example_output --reference /home/Desktop/references/viral_genome_reference.fasta --primer-scheme /home/Desktop/my_primers.bed --adapters /home/Desktop/trimmomatic_adapters/adapters.fa --threads 2 --threads-total 4

The output directory will contain subdirectories with QC data and reports, logs and all files related to the assembly pipeline, including consensus sequences and intermediate files (mapping files, coverage reports).

Tip: Users are strongly encouraged to always use absolute paths. This minimizes the chances for mistakes and enforces the usage of all correct files. 

### Taxonomic assignment

The **metagenomics pipeline** (`viralunity meta`) takes raw sequencing reads to taxonomic classifications and visualizations. You can use **Kraken2 only**, **Diamond only**, or **both**. Each tool can be run on **reads** and, when assembly is enabled, on **contigs**.

#### Sample sheet format

The `--sample-sheet` argument must point to a CSV file (no header row) with sample IDs and file paths.

- **Illumina:** three columns — sample name, R1 path, R2 path.
- **Nanopore:** two columns — sample name, path to single FASTQ/FASTA.

Example Illumina:

    sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
    sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz

Example Nanopore:

    barcode01,/path/to/barcode01.fastq
    barcode02,/path/to/barcode02.fastq.gz

You can create an Illumina sample sheet with `viralunity_create_samplesheet.py` (see `--help`). Prefer absolute paths.

#### Pipeline overview (Illumina)

1. **Quality control** — [fastp](https://github.com/OpenGene/fastp): trimming, quality filtering, adapter handling (auto-detect or custom FASTA), JSON/HTML reports.
2. **Optional dehosting** — If a host reference is provided, [minimap2](https://github.com/lh3/minimap2) aligns reads and unmapped reads are kept for downstream steps.
3. **Merged reads** — Host-filtered (or trimmed) R1 and R2 are merged per sample for read-level classification.
4. **Read classification** — [Kraken2](https://github.com/DerrickWood/kraken2) and/or [DIAMOND](https://github.com/bbuchfink/diamond) blastx on reads (optional each).
5. **Optional de novo assembly** — [MEGAHIT](https://github.com/voutcn/megahit) on host-filtered pairs.
6. **Contig classification** — Kraken2 and/or Diamond on contigs (when assembly is run). For Diamond contigs: viral contigs are extracted, reads are remapped for support, and results are filtered by read count before summarization.
7. **Summaries and filters** — Per-sample taxa tables (family/genus/species), RPM normalization, optional max-RPM bleed filter and negative-control background filter. [Krona](https://github.com/marbl/Krona) plots and [MultiQC](https://github.com/ewels/MultiQC) reports.

Output directories include `qc/` (fastp + FastQC + MultiQC), `host_filtered/`, `metagenomics/taxonomic_assignments/` (per-tool and per-mode results), `denovo_assembly/` and `mapping/viral/` when assembly is enabled.

#### Pipeline overview (Nanopore)

1. **Optional dehosting** — If `--host-reference` is set, minimap2 aligns reads; unmapped reads are kept (no fastp QC).
2. **Read classification** — Kraken2 and/or Diamond on (host-filtered or raw) reads.
3. **Optional de novo assembly** — MEGAHIT on host-filtered reads.
4. **Optional polishing** — [Racon](https://github.com/lbcb-sci/racon) and/or [Medaka](https://github.com/nanoporetech/medaka) on MEGAHIT contigs (Nanopore-only options).
5. **Contig classification** — Kraken2 and/or Diamond on (polished or raw) contigs; Diamond contigs use read support and filtering as in Illumina.
6. **Summaries and filters** — Same as Illumina: taxa tables, RPM, bleed filter, negative-control filter, Krona plots.

Output layout matches Illumina where applicable; polishing outputs go under `denovo_assembly/megahit/{sample}/` and `medaka_work/`.

#### Databases

- **Kraken2** — Required only if using Kraken2. Pre-built indexes are available [here](https://benlangmead.github.io/aws-indexes/k2). Example (viral DB):

      mkdir -p /path/to/kraken2_database/viral
      cd /path/to/kraken2_database/viral
      wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz
      tar -xzvf k2_viral_20240112.tar.gz

- **Krona** — Required when running any classification (Kraken2 and/or Diamond). Update taxonomy as needed:

      mamba activate viralunity
      rm -rf $CONDA_PREFIX/opt/krona/taxonomy
      mkdir -p /path/to/krona/taxonomy
      ln -s /path/to/krona/taxonomy $CONDA_PREFIX/opt/krona/taxonomy
      ktUpdateTaxonomy.sh /path/to/krona/taxonomy

- **Taxdump** — Required when running any classification. Download [NCBI taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz), extract so that `nodes.dmp` and `names.dmp` are in one directory, and pass that path to `--taxdump`.

- **Diamond** — Required only if using Diamond. Provide a protein FASTA (e.g. RefSeq viral) and build the database (pipeline runs `diamond makedb` if needed). You also need an NCBI **assembly summary** file for taxonomy mapping (e.g. `assembly_summary_refseq.txt` from NCBI).

#### Metagenomics pipeline — full option reference

Run `viralunity meta -h` for the complete list. Below is a structured reference of all options.

**Required (all runs)**

| Option | Description |
|--------|-------------|
| `--data-type` | `illumina` or `nanopore`. |
| `--sample-sheet` | Path to CSV with sample IDs and file paths (see Sample sheet format above). |
| `--config-file` | Path to the YAML config file to be generated (and used when running the workflow). |
| `--output` | Base output directory. Results go under `{output}/{run_name}/`. |

**Run identity and resources**

| Option | Default | Description |
|--------|---------|-------------|
| `--run-name` | `undefined` | Name of the run; output is written to `{output}/{run_name}/`. |
| `--threads` | `1` | Threads per job. |
| `--threads-total` | `1` | Total cores for the workflow (Snakemake parallelization). |
| `--create-config-only` | off | Only generate the config file; do not run the workflow. |

**Classification — Kraken2**

| Option | Default | Description |
|--------|---------|-------------|
| `--kraken2-database` | `NA` | Path to Kraken2 database (required if using Kraken2). |
| `--run-kraken2-reads` | on | Enable Kraken2 on reads. |
| `--no-kraken2-reads` | — | Disable Kraken2 on reads. |
| `--run-kraken2-contigs` | on | Enable Kraken2 on contigs when assembly is run. |
| `--no-kraken2-contigs` | — | Disable Kraken2 on contigs. |

**Classification — Diamond**

| Option | Default | Description |
|--------|---------|-------------|
| `--run-diamond-reads` | off | Enable DIAMOND blastx on reads. |
| `--no-diamond-reads` | — | Disable DIAMOND on reads. |
| `--run-diamond-contigs` | off | Enable DIAMOND on assembled contigs (when `--run-denovo-assembly`). |
| `--no-diamond-contigs` | — | Disable DIAMOND on contigs. |
| `--diamond-database` | `NA` | Protein FASTA for Diamond DB (required if running Diamond). |
| `--assembly-summary` | `NA` | NCBI assembly summary for Diamond taxonomy (required if running Diamond). |
| `--diamond-sensitivity` | `sensitive` | One of: `sensitive`, `mid-sensitive`, `more-sensitive`, `ultra-sensitive`. |
| `--evalue` | `0.001` | Diamond E-value threshold. |

**Classification — shared (Kraken2 and/or Diamond)**

| Option | Default | Description |
|--------|---------|-------------|
| `--krona-database` | `NA` | Path to Krona taxonomy (required when running any classification). |
| `--taxdump` | `NA` | Path to NCBI taxdump dir (`nodes.dmp`, `names.dmp`) for taxonomic summaries. |
| `--remove-human-reads` | off | Remove human reads from Krona/summaries. |
| `--remove-unclassified-reads` | off | Remove unclassified reads from Krona/summaries. |

**QC — Illumina only (fastp)**

| Option | Default | Description |
|--------|---------|-------------|
| `--adapters` | `NA` | Path to adapter FASTA; `NA` = auto-detect. |
| `--minimum-read-length` | `50` | Minimum read length after trimming (fastp `--length_required`). |
| `--trim` | `0` | Bases to trim from 5′ end (used as trim_head if `--trim-head` not set). |
| `--trim-head` | — | Bases to trim from 5′ (overrides `--trim`). |
| `--trim-tail` | — | Bases to trim from 3′. |
| `--cut-front-mean-quality` | `20` | fastp cut_front mean quality threshold. |
| `--cut-tail-mean-quality` | `20` | fastp cut_tail mean quality threshold. |
| `--cut-right-window-size` | `4` | fastp cut_right window size. |
| `--cut-right-mean-quality` | `20` | fastp cut_right mean quality threshold. |

**Dehosting and assembly**

| Option | Default | Description |
|--------|---------|-------------|
| `--host-reference` | `NA` | Host genome FASTA for dehosting (minimap2). Omit or `NA` to skip. |
| `--run-denovo-assembly` | off | Run MEGAHIT and classify contigs. |

**Nanopore-only — polishing**

| Option | Default | Description |
|--------|---------|-------------|
| `--run-polish-racon` | off | Run Racon polishing on MEGAHIT assembly. |
| `--no-polish-racon` | — | Disable Racon (default). |
| `--run-polish-medaka` | off | Run Medaka polishing on assembly. |
| `--no-polish-medaka` | — | Disable Medaka (default). |
| `--medaka-model` | — | Medaka model (e.g. `r941_min_high_g360`). Omit to use Medaka default. |

**Summary filters**

| Option | Default | Description |
|--------|---------|-------------|
| `--bleed-fraction` | `0.005` | Max-RPM bleed filter fraction. |
| `--negative-controls` | (empty) | Comma-separated sample IDs to use as negative controls for background filter. |
| `--negative-p-threshold` | `0.01` | p-value threshold for negative-control filter. |

#### How to run

1. Create a sample sheet (see format above) or use `viralunity_create_samplesheet.py` for Illumina.
2. Activate the environment: `conda activate viralunity`.
3. Run `viralunity meta` with `--data-type illumina` or `--data-type nanopore`, plus required and optional arguments below.

The script generates the config file from your arguments and then runs the Snakemake workflow. Use `--create-config-only` to only generate the config and run Snakemake yourself later (e.g. `snakemake --snakefile viralunity/scripts/metagenomics_<illumina|nanopore>.smk --configfile <config.yml> --cores N`).

#### Run examples

**QC only** (fastp + MultiQC, no classification):

    viralunity meta --data-type illumina --sample-sheet samples.csv --config-file config.yaml --output /path/to/output --no-kraken2-reads

**Default: Kraken2 on reads** (QC + classification + Krona + summaries):

    viralunity meta --data-type illumina --sample-sheet samples.csv --config-file config.yaml --output /path/to/output \
      --kraken2-database /path/to/kraken2_db --krona-database /path/to/krona_taxonomy --taxdump /path/to/taxdump

**With adapters and dehosting:**

    viralunity meta --data-type illumina --sample-sheet samples.csv --config-file config.yaml --output /path/to/output \
      --kraken2-database /path/to/kraken2_db --krona-database /path/to/krona_taxonomy --taxdump /path/to/taxdump \
      --adapters /path/to/adapters.fa --host-reference /path/to/host_genome.fa --threads 4 --threads-total 8

**Diamond on reads only** (no Kraken2):

    viralunity meta --data-type illumina --sample-sheet samples.csv --config-file config.yaml --output /path/to/output \
      --no-kraken2-reads --run-diamond-reads --diamond-database /path/to/proteins.faa --assembly-summary /path/to/assembly_summary.txt \
      --krona-database /path/to/krona_taxonomy --taxdump /path/to/taxdump

**Full option: assembly + Kraken2 + Diamond (reads and contigs):**

    viralunity meta --data-type illumina --sample-sheet samples.csv --config-file config.yaml --output /path/to/output \
      --kraken2-database /path/to/kraken2_db --krona-database /path/to/krona_taxonomy --taxdump /path/to/taxdump \
      --run-diamond-reads --run-diamond-contigs --diamond-database /path/to/proteins.faa --assembly-summary /path/to/assembly_summary.txt \
      --run-denovo-assembly --host-reference /path/to/host.fa --threads 4 --threads-total 8

**Nanopore — full pipeline (dehosting, assembly, Medaka polishing, Kraken2 + Diamond):**

    viralunity meta --data-type nanopore --sample-sheet samplesheet_nano.csv --config-file my_config_nano.yml --output output_nano \
      --run-name undefined --threads 4 --threads-total 4 \
      --kraken2-database /path/to/kraken2_db --krona-database /path/to/krona_taxonomy --taxdump /path/to/taxdump \
      --assembly-summary /path/to/assembly_summary.txt --diamond-database /path/to/proteins.faa \
      --host-reference /path/to/host.fa --run-denovo-assembly \
      --run-kraken2-reads --run-kraken2-contigs --run-diamond-reads --run-diamond-contigs \
      --run-polish-medaka --medaka-model r941_min_high_g360 \
      --bleed-fraction 0.005 --negative-p-threshold 0.01

Omit `--run-polish-medaka` and `--medaka-model` if you do not want polishing. Kraken2 on reads/contigs is on by default; add `--no-kraken2-reads` or `--no-kraken2-contigs` to disable.

#### Run

Create a sample sheet (e.g. with `viralunity_create_samplesheet.py` for Illumina), then activate the environment and launch. Use `--data-type illumina` or `--data-type nanopore` and the options from the tables above.

    conda activate viralunity
    viralunity meta --data-type illumina --sample-sheet /path/to/example.csv --config-file /path/to/example_meta.yml \
      --output /path/to/example_output_meta --run-name example_run_meta \
      --kraken2-database /path/to/kraken2_database/viral/ --krona-database /path/to/krona/taxonomy/ --taxdump /path/to/taxdump \
      --threads 2 --threads-total 4

Output includes `qc/` (fastp reports, FastQC, MultiQC — Illumina only), `host_filtered/` (if dehosting is used), `metagenomics/taxonomic_assignments/` (Kraken2 and/or Diamond results, Krona plots, taxa summaries, RPM and filter tables when applicable), and optionally `denovo_assembly/`, `mapping/viral/`, and for Nanopore with polishing `medaka_work/`. 

### Running both pipelines

Users interested in performing both taxonomic assignment followed by consensus sequence inferences are encouraged to store results for both workflows in the same output directory. This allows snakemake to avoid repeating steps present in both pipelines, mainly data quality control.

Currently, the easiest way to do this is by regularly running the metagenomics pipeline (as above) and manually editing the config file for the consensus pipeline. The later config file can be generated using the script viralunity_consensus.py with the --create-config-only, which will prevent the workflow from running automatically. This file should be manually edited, removing samples uninteresting for assembly purpuses and setting the output path to the same one used in the metagenomics workflow. 

Once manual edition is finished, one can execute the workflow directly with snakemake:

    snakemake --snakefile  viralunity/scripts/consensus_<datatype>.smk --configfile /home/Desktop/config_consensus.edited.yml --cores all -p

## Notes

### Nanopore data

The metagenomics pipeline supports Nanopore data with `--data-type nanopore`. There is no fastp-based QC step; you can use optional dehosting, MEGAHIT assembly, and optional Racon/Medaka polishing. See **Pipeline overview (Nanopore)** and **Nanopore-only — polishing** in the option reference above, and the Nanopore run example. Sample sheets for Nanopore have two columns (sample ID, path to one FASTQ/FASTA per sample). 

 
### Segmented viruses

Even though the pipeline has been originally designed to handle non-segmented viruses, it can be naively used to assemble segmented genomes. One just needs to specify one genomic segment as reference at a time. This will automatically create output directory for each segment, which can be analyzed in downstream workflows. This solution is far from optimal, and better arrangements will be available on following versions.

## Tests

To execute unit tests run
```sh
make test
```

## Citation

A scientific publication fully describing this pipeline is being prepared. Meanwhile, feel free to cite this GitHub repo. Primary references for used dependencies should also be cited:

<a href="https://github.com/OpenGene/fastp">fastp</a>: Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018; 34(17):i884–i890.

<a href="http://www.usadellab.org/cms/?page=trimmomatic">Trimmomatic</a>: Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014; 30(15):2114–212.

<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">Fastqc</a>

<a href="https://github.com/ewels/MultiQC">Multiqc</a>: Ewels P, Magnusson M, Lundin S, et al. MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016; 32(19):3047–8.

<a href="https://github.com/lh3/minimap2">Minimap2</a>: Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018; 34:3094-3100.

<a href="https://github.com/samtools/samtools">Samtools</a>: Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009; 25(16):2078–9.

<a href="https://github.com/andersen-lab/ivar">iVar</a>: Grubaugh ND, Gangavarapu K, Quick J, et al. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biol 20, 8 (2019)

<a href="https://github.com/arq5x/bedtools2">BEDtools</a>: Quinlan AR, Hall IM. BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics. 2010; 26(6):841–2.

<a href="https://doi.org/10.1186/s13059-019-1891-0">Kraken2</a>: Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019).

<a href="https://github.com/bbuchfink/diamond">DIAMOND</a>: Buchfink B, Reuter K, Drost H-G. Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods 18, 366–368 (2021).

<a href="https://github.com/voutcn/megahit">MEGAHIT</a>: Li D, Liu C-M, Luo R, Sadakane K, Lam T-W. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics. 2015; 31(10):1674–1676.

<a href="https://doi.org/10.1186/1471-2105-12-385">Krona</a>: Ondov BD, Bergman NH, Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics 12, 385 (2011). 

For visualization of fasta and bam files, we recommend <a href="https://ormbunkar.se/aliview/">Aliview</a> and <a href="https://ics.hutton.ac.uk/tablet/">Tablet</a>, respectively.
