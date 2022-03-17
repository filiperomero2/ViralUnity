# ViralUnity

ViralUnity is yet another simple tool to perform data quality control and viral genome reference assembly from Illumina paired end reads. It is designed to be effective and easy to use.

## Installation

ViralUnity is distributed as single bash script. Nevertheless, it relies on diverse dependencies and all of them have been conveniently documented on a conda environment yaml file. From this file one can easily install required softwares and run the pipeline.

To enable ViralUnity, clone the repo and create the environment:

    $ cd
    $ git clone https://github.com/filiperomero2/ViralUnity.git
    $ cd ~/ViralUnity/ENV/
    $ conda env create --file ViralUnity.yml


## Usage

### Arguments and data organization requirements

ViralUnity was designed to be as simple as possible, with the objective of making users go from raw reads to processed consensus genome sequences for entire batches of samples with a single command line. The script may take six arguments, being only the first four strictly required:

    1 - The absolute path for samples root directory;
    2 - The absolute path for output directory;
    3 - The absolute path for a reference genome in fasta format;
    4 - The absolute path for trimmomatic adapters fasta file;
    5 - the minimum sequencing depth necessary to incorporate a base into the consensus sequence (default: 100);
    6 - the number of threads available for processing (default: 1).

While arguments 2 to 6 are self-explanatory, argument 1 may demand clarification. To be able to analyze entire sequencing runs from a single command line, ViralUnity needs data to be stored in a well specified structure of directories. The path in argument 1 refers to a directory that harbors samples' directories, each containing two fastq files (R1 and R2 reads), like in the example bellow:

    (base) fmoreira@DESKTOP:~/Desktop/pilot/DATA$ tree
    .
    ├── SEQ_176_L001_ds.5c9c29d7d8824301ae5a4ddee48ee059
    │   ├── SEQ-176_S74_L001_R1_001.fastq.gz
    │   └── SEQ-176_S74_L001_R2_001.fastq.gz
    └── SEQ_180_L001_ds.b269a1119df645ad8c2876cd312620a3
        ├── SEQ-180_S35_L001_R1_001.fastq.gz
        └── SEQ-180_S35_L001_R2_001.fastq.gz

    2 directories, 4 files

Importantly, in case this directory organization is not strictly available, the pipeline will fail. Any directory outside these specifications on the root directory will halt execution. Each sample directory is expected to only have two (R1 and R2) fastq files (with .fastq or .fastq.gz extensions). Sample names should not include the underline character (_). 

If the structure is correct, the script iterates over these directories, performing QC (trimmomatic), read mapping (Bowtie2), variant calling (bcftools) and consensus sequence inferences (bcftoools/bedtools). QC reports are generated with fastQC and multiQC. The results from all these analysis are stored in the specified output directory. 

### Run

To run the pipeline, just activate the conda environment and launch the analysis:

    $ conda activate ViralUnity
    $ ~/ViralUnity/SCRIPT/ViralUnity.sh --LIBDIR ~/LIBRARIES/RUN_1/ --OUTDIR ~/ANALYSIS/RUN1/ --REF ~/REFERENCE_GENOMES/reference.fasta --ADAPTERS ~/trimmomatic/adapter.fa

One may also specify custom depth thresholds and number of processing threads:

    $ ~/ViralUnity/SCRIPT/ViralUnity.sh --LIBDIR ~/LIBRARIES/RUN_1/ --OUTDIR ~/ANALYSIS/RUN1/ --REF ~/REFERENCE_GENOMES/reference.fasta --ADAPTERS ~/trimmomatic/adapter.fa --MINCOV 200 --THREADS 6

The output directory contains 2 report files, one with assembly statistics and other with a timestamp for each sample processing. In addition, four directories are also created, comprehending QC reports for raw and filtered data, mapping and variants associated files and consensus sequences. 

### Note on segmented viruses

Even though the pipeline has been originally designed to handle non-segmented viruses, it can be naively used to assemble segmented genomes. One just needs to specify one genomic segment as reference at a time. This will automatically create output directory for each segment, which can be analyzed in downstream workflows.

### Note on primer sequences removal

The removal of primer associated SNPs is a mandatory step in processing sequences generated from targetted sequencing approaches. While this pipeline does not use any tool that specifically query for primer sequences, it does crop out the initial 30 bp of all reads with trimmomatic. As primer sequences hardly extends more than 30 bp, this step avoids the introduction of artefactual SNPs. Users analyzing data generated under alternative protocols (e.g., metagenomics), may want to change this behavior by removing the HEADCROP:30 trimmomatic argument in the script.


## Citation

If you use this pipeline, please cite this github repo and also: 

<a href="http://www.usadellab.org/cms/?page=trimmomatic">Trimmomatic</a>: Anthony M. Bolger, Marc Lohse, Bjoern Usadel, Trimmomatic: a flexible trimmer for Illumina sequence data, Bioinformatics, Volume 30, Issue 15, 1 August 2014, Pages 2114–212.

<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">Fastqc</a>

<a href="https://github.com/ewels/MultiQC">Multiqc</a>: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047–8.

<a href="https://github.com/BenLangmead/bowtie2">Bowtie2</a>: Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012;9(4):357–9.

<a href="https://github.com/samtools/samtools">Samtools</a>: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078–9.

<a href="https://github.com/vcftools/vcftools">VCFtools</a>: Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, et al. The variant call format and VCFtools. Bioinformatics. 2011;27(15):2156–8.

<a href="https://github.com/arq5x/bedtools2">BEDtools</a>: Quinlan AR, Hall IM. BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841–2.

For visualization of fasta and bam files, we recommend <a href="https://ormbunkar.se/aliview/">Aliview</a> and <a href="https://ics.hutton.ac.uk/tablet/">Tablet</a>, respectively.
