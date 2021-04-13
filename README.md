# ViralUnit

ViralUnit is a tool to perform quality control and viral genome reference assembly from Illumina paired end reads. While other pipelines have been developed for this purpose, we believe ViralUnit has a unique combination of efficiency, scalability and ease of use. 

## Instalation

ViralUnit is distributed as single bash script. Nevertheless, it relies on diverse dependencies to do its job. All these softwares have been conveniently documented on a single yaml file, suitable for construction of a conda environment. To do this, clone or download this repo and simply:

    $ cd; git clone https://github.com/filiperomero2/ViralUnit.git
    $ cd ~/ViralUnit/ENV/
    $ conda create env --file ViralUnit.yml

These commands will automatically create an environment with all necessary dependencies and the pipeline will be enabled.

## Usage

### Arguments and data organization requirements

ViralUnit was designed to be as simple as possible, with the objective of making users go from raw reads to processed consensus genome sequences for entire batches of samples with a single command line. The script may take four positional arguments, being only the first two strictly required:

    1 - The path for samples' root directory;
    2 - The path for a reference genome in fasta format;
    3 - the minimum sequencing depth necessary to incorporate a base into the consensus sequence (default: 100);
    4 - the number of threads available for processing (default: 1).

While arguments 2, 3 and 4 are self-explanatory, argument 1 may demand clarification. To be able to analyze entire sequencing runs from a single command line, ViralUnit needs data to be stored in a well specified structure of directories. The path in argument 1 refers to a directory that harbors samples' directories, each containing two fastq files (R1 and R2 reads), like in the example bellow:

    (base) fmoreira@DESKTOP-IIHQFUE:~/Desktop/pilot/DATA$ tree
    .
    ├── SEQ_176_L001_ds.5c9c29d7d8824301ae5a4ddee48ee059
    │   ├── SEQ-176_S74_L001_R1_001.fastq.gz
    │   └── SEQ-176_S74_L001_R2_001.fastq.gz
    └── SEQ_180_L001_ds.b269a1119df645ad8c2876cd312620a3
        ├── SEQ-180_S35_L001_R1_001.fastq.gz
        └── SEQ-180_S35_L001_R2_001.fastq.gz

    2 directories, 4 files

The script then iterates over these directories, performing QC, read mapping, variant calling and consensus sequence inferences. The results from all these analysis is stored in a directory with the suffix RESULTS/, created within the path specified by argument 1. 

### Run

To effectively run the pipeline, just activate the conda environment and use a simple command line, examplified as follows:

    $ conda activate ViralUnit
    $ ./ViralUnit.sh ~/LIBRARIES/RUN_1/ ~/REFERENCE_GENOMES/reference.fasta 200 6

These lines would make all dependencies available, start the analysis for all samples stored in the first path, using as reference the fasta file specified in the second. Additionaly, consensus genome sequences sites with less then 200x depth would be masked and the analysis would use 6 threads. 

### Note on segmented viruses

Even though the pipeline has been originally designed to handle non-segmented viruses, it can be naively used to assemble segmented genomes. One just needs to specify one genomic segment as reference at a time. This will automatically create a RESULTS/ directory for each segment, which can be jointly analyzed in downstream workflows.

## Citation

If you use this pipeline, please cite: 

<a href="https://github.com/OpenGene/fastp">Fastp</a>: Chen S, Zhou Y, Chen Y, Gu J. Fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17): 884–90.

<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">Fastqc</a>

<a href="https://github.com/ewels/MultiQC">Multiqc</a>: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047–8.

<a href="https://github.com/BenLangmead/bowtie2">Bowtie2</a>: Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012;9(4):357–9.

<a href="https://github.com/samtools/samtools">Samtools</a>: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078–9.

<a href="https://github.com/vcftools/vcftools">VCFtools</a>: Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, et al. The variant call format and VCFtools. Bioinformatics. 2011;27(15):2156–8.

<a href="https://github.com/arq5x/bedtools2">BEDtools</a>: Quinlan AR, Hall IM. BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841–2.
