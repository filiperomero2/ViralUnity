# ViralUnity

ViralUnity is a simple tool to perform analysis of viral high-throughput sequencing data. It comprises a collection of python scripts and snakemake workflows, including several software dependecies, to perform data quality control, taxonomic assignments and reference genome assembly.

ViralUnity runs on *nix systems and is able to process entire sequencing runs in minimal time on a regular computer.

## Installation

ViralUnity is python package that launches snakemake workflows. All dependencies are documented in the conda environment file `environment.yml`. We recommend using <a href="https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html">mamba</a> for installation. 

To enable ViralUnity, clone the repo and create the environment:

    git clone https://github.com/filiperomero2/ViralUnity.git
    cd viralunity/
    conda create -n viralunity python=3.11
    conda activate viralunity
    conda env update -n viralunity --file environment.yml
    pip install .

Tip: On macOS (with chips m1 or later), one may need to further configure the environment before installing dependencies:

    conda create -n viralunity
    conda activate viralunity
    conda config --env --set subdir osx-64
    conda env update -n viralunity --file environment.yml
    pip install .

It is a good practice to verify that all required dependencies (listed in `environment.yml`) are available before running the pipeline. On certain systems it may be necessary to manually set conda channel priority.

## Usage

ViralUnity comprehends two main pipelines, embodied in separate scripts and snakemake workflows, for performing reference genome assembly (viralunity_consensus.py) and reads taxonomic assignment (viralunity_meta.py). 

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

The output directory will contain subdirectories with QC data and reports, logs and all files related to the asssembly pipeline, including consensus sequences and intermediate files (mapping files, coverage reports). 

Tip: Users are strongly encouraged to always use absolute paths. This minimizes the chances for mistakes and enforces the usage of all correct files. 

### Taxonomic assignment

The viralunity metagenomics pipeline allows users to go from raw sequencing reads to metagenomics classifications and visualizations, characterizing viral diversity. You can choose **Kraken2 only**, **Diamond only**, or **both** for classification; each tool can be run on reads and (optionally) on assembled contigs.

#### Pipeline overview (Illumina)

1. **Quality control** — [fastp](https://github.com/OpenGene/fastp): trimming, quality filtering, adapter handling (auto-detect or custom FASTA), JSON/HTML reports.
2. **Optional dehosting** — If a host reference is provided, [minimap2](https://github.com/lh3/minimap2) aligns reads and unmapped reads are kept for downstream steps.
3. **Merged reads** — Host-filtered (or trimmed) R1 and R2 are merged per sample for read-level classification.
4. **Read classification** — [Kraken2](https://github.com/DerrickWood/kraken2) and/or [DIAMOND](https://github.com/bbuchfink/diamond) blastx on reads (optional each).
5. **Optional de novo assembly** — [MEGAHIT](https://github.com/voutcn/megahit) on host-filtered pairs.
6. **Contig classification** — Kraken2 and/or Diamond on contigs (when assembly is run). For Diamond contigs: viral contigs are extracted, reads are remapped for support, and results are filtered by read count before summarization.
7. **Summaries and filters** — Per-sample taxa tables (family/genus/species), RPM normalization, optional max-RPM bleed filter and negative-control background filter. [Krona](https://github.com/marbl/Krona) plots and [MultiQC](https://github.com/ewels/MultiQC) reports.

Output directories include `qc/` (fastp + FastQC + MultiQC), `host_filtered/`, `metagenomics/taxonomic_assignments/` (per-tool and per-mode results), `denovo_assembly/` and `mapping/viral/` when assembly is enabled.

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

#### Metagenomics arguments

**Required:** `--data-type`, `--sample-sheet`, `--config-file`, `--output`

**Classification (choose one or both):** `--kraken2-database`, `--krona-database`, `--taxdump` (all three required when running any classification); `--no-kraken2-reads`, `--no-kraken2-contigs`; `--run-diamond-reads`, `--run-diamond-contigs`, `--diamond-database`, `--assembly-summary`; `--diamond-sensitivity`, `--evalue`

**QC (fastp, Illumina):** `--adapters` (optional, default NA = auto-detect), `--minimum-read-length`, `--trim`, `--trim-head`, `--trim-tail`, `--cut-front-mean-quality`, `--cut-tail-mean-quality`, `--cut-right-window-size`, `--cut-right-mean-quality`

**Dehosting and assembly:** `--host-reference`, `--run-denovo-assembly`

**Filters and other:** `--remove-human-reads`, `--remove-unclassified-reads`, `--bleed-fraction`, `--negative-controls`, `--negative-p-threshold`, `--run-name`, `--create-config-only`, `--threads`, `--threads-total`

Run `viralunity meta -h` for full help.

#### Minimal run examples

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

#### Run

Create a sample sheet (e.g. with `viralunity_create_samplesheet.py`), then activate the environment and launch:

    conda activate viralunity
    viralunity meta --data-type illumina --sample-sheet /path/to/example.csv --config-file /path/to/example_meta.yml \
      --output /path/to/example_output_meta --run-name example_run_meta \
      --kraken2-database /path/to/kraken2_database/viral/ --krona-database /path/to/krona/taxonomy/ --taxdump /path/to/taxdump \
      --threads 2 --threads-total 4

Output includes `qc/` (fastp reports, FastQC, MultiQC), `host_filtered/` (if dehosting is used), `metagenomics/taxonomic_assignments/` (Kraken2 and/or Diamond results, Krona plots, taxa summaries, RPM and filter tables when applicable), and optionally `denovo_assembly/` and `mapping/viral/`. 

### Running both pipelines

Users interested in performing both taxonomic assignment followed by consensus sequence inferences are encouraged to store results for both workflows in the same output directory. This allows snakemake to avoid repeating steps present in both pipelines, mainly data quality control.

Currently, the easiest way to do this is by regularly running the metagenomics pipeline (as above) and manually editing the config file for the consensus pipeline. The later config file can be generated using the script viralunity_consensus.py with the --create-config-only, which will prevent the workflow from running automatically. This file should be manually edited, removing samples uninteresting for assembly purpuses and setting the output path to the same one used in the metagenomics workflow. 

Once manual edition is finished, one can execute the workflow directly with snakemake:

    snakemake --snakefile  viralunity/scripts/consensus_<datatype>.smk --configfile /home/Desktop/config_consensus.edited.yml --cores all -p

## Notes

### Nanopore data

Experimental pipelines for nanopore sequencing data were recently added. It is still a work in progress, and there's a lot of room for improvements. Unlike the Illumina pipeline, it does not include data quality control steps. In case you want to use it, just use the flag '--data-type nanopore'. 

 
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
