# ViralUnity

ViralUnity is yet another simple tool to perform data quality control and viral genome reference assembly from Illumina paired-end reads. It is designed to be effective and easy to use. It runs on *nix systems and is able to process complete Illumina MiSeq runs in a couple of hours on a regular laptop.

## Installation

ViralUnity is distributed as single bash script. Nevertheless, it relies on diverse dependencies and all of them have been conveniently documented on a conda environment YAML file. From this file one can easily install required softwares and run the pipeline.

To enable ViralUnity, clone the repo and create the environment:

    $ cd
    $ git clone https://github.com/filiperomero2/ViralUnity.git
    $ cd ~/ViralUnity/ENV/
    $ conda env create --file ViralUnity.yml

In case of compatibility issues between OS and versions specified in the yaml file, one may create the environment and download dependencies manually from the bioconda channel:

    $ conda create --name ViralUnity
    $ conda activate ViralUnity
    $ conda install -c bioconda fastqc multiqc trimmomatic minimap2 samtools bcftools bedtools -y

Keep in mind this pipeline has been developed and tested on the following dependencies versions:

* fastqc v0.11.9
* multiqc v1.9
* trimmomatic v0.39
* minimap2 v2.17
* samtools v1.11
* bcftools v1.11
* bedtools v2.30.0

## Usage

### Arguments and data organization requirements

ViralUnity was designed to be as simple as possible, with the objective of making users go from raw reads to processed consensus genome sequences for entire batches of samples with a single command line. The script may take eight arguments, being only the first four strictly required:

    --LIBDIR - The absolute path for samples root directory;
    --OUTDIR - The absolute path for output directory;
    --REF - The absolute path for a reference genome in fasta format;
    --ADAPTERS - The absolute path for trimmomatic adapters fasta file;
    --MINCOV - the minimum sequencing depth necessary to incorporate a base into the consensus sequence (default: 100);
    --MINLEN - Minimum read length (Optional; default = 50);
    --HEADCROP - Number of bases to trim from the start of the read, useful for primer sequences removal (Optional; default = 30);
    --THREADS - the number of threads available for processing (default: 1).

The first argument may demand clarification. To be able to analyze entire sequencing runs from a single command line, ViralUnity needs data to be stored in a well specified structure of directories. The path in argument 1 refers to a directory that harbors samples' directories, each containing two fastq files (R1 and R2 reads), like in the example bellow:

    (base) username@DESKTOP:~/Desktop/pilot/DATA$ tree
    .
    ├── SEQ_176_L001_ds.5c9c29d7d8824301ae5a4ddee48ee059
    │   ├── SEQ-176_S74_L001_R1_001.fastq.gz
    │   └── SEQ-176_S74_L001_R2_001.fastq.gz
    └── SEQ_180_L001_ds.b269a1119df645ad8c2876cd312620a3
        ├── SEQ-180_S35_L001_R1_001.fastq.gz
        └── SEQ-180_S35_L001_R2_001.fastq.gz

    2 directories, 4 files

Importantly, in case this directory organization is not strictly available, the pipeline will fail. Any directory outside these specifications on the root directory will halt execution. Each sample directory is expected to only have two (R1 and R2) fastq files (with .fastq or .fastq.gz extensions). Sample names should not include the underline character (_). 

If the structure is correct, the script iterates over these directories, performing QC (trimmomatic), read mapping (minimap2), variant calling (bcftools) and consensus sequence inferences (bcftoools/bedtools). QC reports are generated with fastQC and multiQC. The results from all these analysis are stored in the specified output directory. 

### Run

To run the pipeline, just activate the conda environment and launch the analysis:

    $ conda activate ViralUnity
    $ ~/ViralUnity/SCRIPT/ViralUnity.sh --LIBDIR ~/LIBRARIES/RUN_1/ --OUTDIR ~/ANALYSIS/RUN1/ --REF ~/REFERENCE_GENOMES/reference.fasta --ADAPTERS ~/trimmomatic/adapter.fa

One may also specify other parameters:

    $ ~/ViralUnity/SCRIPT/ViralUnity.sh --LIBDIR ~/LIBRARIES/RUN_1/ --OUTDIR ~/ANALYSIS/RUN1/ --REF ~/REFERENCE_GENOMES/reference.fasta --ADAPTERS ~/trimmomatic/adapter.fa --MINCOV 200 --MINLEN 30 --HEADCROP 20 --THREADS 6

The output directory contains 2 report files, one with assembly statistics and other with a timestamp for each sample processing. In addition, three directories are also created, comprehending QC reports for raw and filtered data, mapping and variants associated files and consensus sequences. 

## Notes

### Note on segmented viruses

Even though the pipeline has been originally designed to handle non-segmented viruses, it can be naively used to assemble segmented genomes. One just needs to specify one genomic segment as reference at a time. This will automatically create output directory for each segment, which can be analyzed in downstream workflows.

### Note on primer sequences removal

The removal of primer associated SNPs is a mandatory step in processing sequences generated from targetted sequencing approaches. While this pipeline does not use any tool that specifically query for primer sequences, its default mode trim out the initial 30 bp of all reads with trimmomatic. As primer sequences hardly extends for more than 30 bp, this step avoids the introduction of artefactual SNPs. Users analyzing data generated under alternative protocols (e.g., metagenomics) may want to change this behavior by setting the value of --HEADCROP to 0 (zero).


## Citation

If you use this pipeline, please cite this github repo and also: 

<a href="http://www.usadellab.org/cms/?page=trimmomatic">Trimmomatic</a>: Anthony M. Bolger, Marc Lohse, Bjoern Usadel. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014; 30(15):2114–212.

<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">Fastqc</a>

<a href="https://github.com/ewels/MultiQC">Multiqc</a>: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016; 32(19):3047–8.

<a href="https://github.com/lh3/minimap2">Minimap2</a>: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100.

<a href="https://github.com/samtools/samtools">Samtools</a>: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009; 25(16):2078–9.

<a href="https://samtools.github.io/bcftools/">BCFtools</a>: Heng Li. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011; 27(21):2987–2993. 

<a href="https://github.com/arq5x/bedtools2">BEDtools</a>: Quinlan AR, Hall IM. BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics. 2010; 26(6):841–2.

For visualization of fasta and bam files, we recommend <a href="https://ormbunkar.se/aliview/">Aliview</a> and <a href="https://ics.hutton.ac.uk/tablet/">Tablet</a>, respectively.
