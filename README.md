# ViralUnity

ViralUnity is a simple tool to perform analysis of viral high-throughput sequencing data, especifically Illumina paired-end reads. It comprises a collection of python scripts and snakemake workflows, including several software dependecies, to perform data quality control, taxonomic assignments and reference genome assembly.

ViralUnity runs on *nix systems and is able to process entire sequencing runs in minimal time on a regular computer.

## Installation

ViralUnity is a collection of python scripts and snakemake workflows. All dependencies have been conveniently documented on a conda environment file (vu_dependencies.yml). From this file, one can easily install required softwares to run the pipeline. We recommend using <a href="https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html">mamba</a> for this purpose. 

To enable ViralUnity, clone the repo and create the environment:

    git clone https://github.com/filiperomero2/ViralUnity.git
    cd viralunity/
    mamba env create -n viralunity
    mamba env update -n viralunity --file envs/vu_dependencies.yml

Tip: On macOS (with chips m1 or later), one may need to further configure the environment before installing dependencies:

    mamba create -n viralunity
    mamba activate viralunity
    conda config --env --set subdir osx-64
    mamba env update -n viralunity --file envs/vu_dependencies.yml

It is a good practice to verify whether all required dependencies (available on envs/vu_dependecies.yml) are available before running the pipeline. In case any installation fail, try re-installing separetely on the command line with mamba. On certain systems, it might be necessary to manually set conda channels priority.

## Usage

ViralUnity comprehends two main pipelines, embodied in separate scripts and snakemake workflows, for respectivelly performing reference genome assembly (viralunity_consensus.py) and reads taxonomic assignment (viralunity_meta.py). 

### Reference viral genome assembly

The viralunity consensus sequence pipeline was designed to be as simple as possible, with the objective of making users go from raw reads to processed consensus genome sequences with a single command line. The core script is called viralunity_consensus.py and accepts several arguments:

    --sample-sheet          Complete path for a csv file with samples data paths and metadata
    --config-file           Complete path for input (viralunity config) file to be created.
    --output                Complete path for output directory to be created by viralunity.
    --run-name              Name for the sequencing run (optional).
    --reference             Complete path for the reference genome in fasta format
    --adapters              Complete path for Illumina adapters sequences in fasta format
    --minimum-read-length   Minimum read length threshold (Default = 50)
    --minimum-coverage      Minimum sequencing coverage for including base in consensus sequence (Default = 20)
    --trim                  Number of bases to trim from the 5' end of reads (Default = 0)
    --create-config-only    Only create config file, not running the workflow (boolean)
    --threads               Number of available threads for individual tasks (Default = 1)
    --threads-total         Number of available threads for the entire workflow (Default = 1)
    -h, --help              Show this help message and exit

The first option receives as argument the path for a samplesheet (csv file) that contains sample names and fastq file paths. This file can be created automatically with the script create_viralunity_samplesheet.py (see its --help option). 

If the specified files paths are correct, the script will generate config file containing all information required to execute the snakemake consensus workflow, allowing optimal parallelization schemes. Beyond creating the config file, the script will execute the workflow, performing data quality control (trimmomatic), read mapping (minimap2), and consensus sequence inference (samtools). QC reports are generated with fastQC and multiQC. The results from all these analysis are stored in the specified output directory. 

#### Run

To run the pipeline, go to the repository directory,activate the conda environment and create a samplesheet file:

    python scripts/create_viralunity_samplesheet.py --input /home/Desktop/test_data/ --output /home/Desktop/example.csv

Check the contents of the samplesheet file. If the correct paths for fastq files are specified, launch the analysis:

    conda activate viralunity
    cd viralunity/
    python scripts/viralunity_consensus.py --sample-sheet /home/Desktop/example.csv --config-file /home/Desktop/example.yml --run-name example_run --output /home/Desktop/example_output --reference /home/Desktop/references/viral_genome_reference.fasta --adapters /home/Desktop/trimmomatic_adapters/adapters.fa --trim 30 --threads 2 --threads-total 4

The output directory will contain subdirectories with QC data and reports, logs and all files related to the asssembly pipeline, including consensus sequences and intermediate files (mapping and variant calls files, coverage reports). 

Tip: Users are strongly encouraged to always use absolute paths. This minimizes the chances for mistakes and enforces the usage of all correct files. 

### Taxonomic assignment

The viralunity metagenomics pipeline is a more recent development, and is also designed to be as simple as possible. It allows users to go from raw sequencing reads to metagenomics classifications and visualizations, characterizing viral diversity captured through sequencing. It is by far not the most complex (or exhaustive) metagenomics workflow, but it is purpose is to provide a characterization of underlying viral diversity (including possible pathogens) as quickly as possible. 

An additional complexity of running this pipeline revolves around downloading and correctly setting databases. 

Kraken2 databases are available <a href="https://benlangmead.github.io/aws-indexes/k2">here</a>. For instance, to install the viral kraken2 database, one can execute:

    mkdir -p /home/Desktop/kraken2_database/viral
    cd /home/Desktop/kraken2_database/viral
    # check updates whenever installing
    wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz
    tar -xzvf k2_viral_20240112.tar.gz
    cd

For krona, a standard outdated taxonomic database is installed in the conda environment directory. To make updates easier, we encourage its exclusion, followed by the creation of a symbolically linked database elsewhere. To do so, execute:
  
    # tip: exact paths are displayed on terminal at the moment of installation with mamba
    mamba activate viralunity
    rm -rf /home/path/to/conda/envs/viralunity/opt/krona/taxonomy
    mkdir -p /home/Desktop/krona/taxonomy
    ln -s ~/home/Desktop/krona/taxonomy /home/path/to/conda/envs/viralunity/opt/krona/taxonomy
    ktUpdateTaxonomy.sh /home/Desktop/krona/taxonomy

The core script is called viralunity_meta.py and accepts several arguments:

     --sample-sheet                Complete path for a csv file with samples data paths and metadata
     --config-file                 Complete path for input (viralunity config) file to be created.
     --output                      Complete path for output directory to be created by viralunity.
     --run-name                    Name for the sequencing run (optional).
     --kraken2-database            Complete path for the kraken2 database directory
     --krona-database              Complete path for the krona taxonomic database
     --adapters                    Complete path for Illumina adapters sequences in fasta format
     --minimum-read-length         Minimum read length threshold (Default = 50)
     --trim                        Number of bases to trim from the 5' end of reads (Default = 0)
     --remove-human-reads          Remove human reads from krona plot (boolean)
     --remove-unclassified-reads   Remove unclassified reads from krona plot (boolean)
     --create-config-only          Only create config file, not running the workflow (boolean)
     --threads                     Number of available threads for individual tasks (Default = 1)
     --threads-total               Number of available threads for the entire workflow (Default = 1)
     -h, --help                    Show this help message and exit

If all paths are correctly set, the script will generate a config file and run the snakemake metagenomics workflow. Briefly, the pipeline will perform data quality control (trimmomatic), taxonomic assignment (kraken2) and generate visualizations (krona).

#### Run

As performed for the previous workflow, to run the metagenomics pipeline one needs to create a samplesheet file and activate the conda environment. To launch the analysis:

    conda activate viralunity    
    cd viralunity/
    python scripts/viralunity_meta.py --sample-sheet /home/Desktop/example.csv --config-file /home/Desktop/example_meta.yml --run-name example_run_meta --kraken2-database /home/Desktop/kraken2_database/viral/ --krona-database /home/Desktop/krona/taxonomy/ --adapters /home/Desktop/trimmomatic_adapters/adapters.fa --threads 2 --threads-total 4 --output /home/Desktop/example_output_meta

Three directories are created, respectivelly containing logs, quality control data and report, and metagenomics results (kraken2 reports and krona interactive plots).The later directory will also contain one summary file containing information of the viral diversity captured at the family level, as calculated by kraken2. 

### Running both pipelines

Users interested in performing both taxonomic assignment followed by consensus sequence inferences are encouraged to store results for both workflows in the same output directory. This allows snakemake to avoid repeating steps present in both pipelines, mainly data quality control.

Currently, the easiest way to do this is by regularly running the metagenomics pipeline (as above) and manually editing the config file for the consensus pipeline. The later config file can be generated using the script viralunity_consensus.py with the --create-config-only, which will prevent the workflow from running automatically. This file should be manually edited, removing samples uninteresting for assembly purpuses and setting the output path to the same one used in the metagenomics workflow. 

Once manual edition is finished, one can execute the workflow directly with snakemake:

    snakemake --snakefile  scripts/consensus.smk --configfile /home/Desktop/config_consensus.edited.yml --cores all -p

## Notes

### Legacy

Previous versions of the pipeline were released as bash scripts. Despite the ease of use, these versions were remarkbly less flexible and slower than the current snakemake workflow implementation. Notwithstanding, the core legacy script should still work and is available in this repository (legacy/viralunity.sh). 
 
### Segmented viruses

Even though the pipeline has been originally designed to handle non-segmented viruses, it can be naively used to assemble segmented genomes. One just needs to specify one genomic segment as reference at a time. This will automatically create output directory for each segment, which can be analyzed in downstream workflows. This solution is far from optimal, and better arrangements will be available on following versions.

### Primer sequences removal

The removal of primer associated SNPs is a mandatory step in processing sequences generated from targeted sequencing approaches. While this pipeline does not use any tool that specifically query for primer sequences, it offers a general trimming functionality with trimmomatic. For analyzing data generate under an amplicon sequencing schemes, we suggest '--trim 30' as a sensible choice.

### Adding scripts to path

Users might want to use scripts without fully specifying their path. This can be acomplished by adding the scripts directory to the system path. This may vary among systems, but for linux (ubuntu-20.04), the following should work: 

    export PATH="/path/to/viralunity/scripts/:$PATH"

Notice that this will work only for a given terminal session. To make it permanent, this line must be included at the end of the .bashrc file. Also, any change in the scripts/ path will cause this to brake. 

## Citation

A scientific publication fully describing this pipeline is being prepared. Meanwhile, feel free to cite this GitHub repo. Primary references for used dependecies should also be cited:

<a href="http://www.usadellab.org/cms/?page=trimmomatic">Trimmomatic</a>: Anthony M. Bolger, Marc Lohse, Bjoern Usadel. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014; 30(15):2114–212.

<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">Fastqc</a>

<a href="https://github.com/ewels/MultiQC">Multiqc</a>: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016; 32(19):3047–8.

<a href="https://github.com/lh3/minimap2">Minimap2</a>: Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018; 34:3094-3100.

<a href="https://github.com/samtools/samtools">Samtools</a>: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009; 25(16):2078–9.

<a href="https://samtools.github.io/bcftools/">BCFtools</a>: Li, H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011; 27(21):2987–2993. 

<a href="https://github.com/arq5x/bedtools2">BEDtools</a>: Quinlan AR, Hall IM. BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics. 2010; 26(6):841–2.

<a href="https://doi.org/10.1186/s13059-019-1891-0">kraken2</a> - Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). 

<a href="https://doi.org/10.1186/1471-2105-12-385">krona</a> - Ondov, B.D., Bergman, N.H. & Phillippy, A.M. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics 12, 385 (2011). 

For visualization of fasta and bam files, we recommend <a href="https://ormbunkar.se/aliview/">Aliview</a> and <a href="https://ics.hutton.ac.uk/tablet/">Tablet</a>, respectively.