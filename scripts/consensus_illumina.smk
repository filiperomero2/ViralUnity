rule all:
    input:
        config['output'] + "assembly/consensus/final_consensus/aln.consensus.fasta"

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample].split()
    return(reads)

rule perform_qc:
    input:
        get_map_input_fastqs,
    output:
        paired_R1 = config['output'] + "qc/data/trim.p.{sample}_R1.fastq.gz",
        unpaired_R1 = config['output'] + "qc/data/trim.u.{sample}_R1.fastq.gz",
        paired_R2 = config['output'] + "qc/data/trim.p.{sample}_R2.fastq.gz",
        unpaired_R2 = config['output'] + "qc/data/trim.u.{sample}_R2.fastq.gz"
    params:
        minimum_length = config["minimum_length"],
        adapters = config["adapters"],
        trim = config["trim"]
    threads: config["threads"]
    log:
        config['output'] + "logs/trimmomatic/{sample}.log"
    benchmark:
        config['output'] + "logs/trimmomatic/{sample}.benchmark.txt"
    shell:
        "trimmomatic PE -quiet -threads {threads} -phred33 {input} "
        "{output.paired_R1} {output.unpaired_R1} "
        "{output.paired_R2} {output.unpaired_R2} "
        "ILLUMINACLIP:{params.adapters}:2:30:10 "
        "LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 "
        "HEADCROP:{params.trim} MINLEN:{params.minimum_length} 2> {log}"

rule generate_qc_report:
    input:
        config['output'] + "qc/data/trim.p.{sample}_R1.fastq.gz",
        config['output'] + "qc/data/trim.p.{sample}_R2.fastq.gz"
    output:
        config['output'] + "qc/reports/trim.p.{sample}_R1_fastqc.html",
        config['output'] + "qc/reports/trim.p.{sample}_R2_fastqc.html",
        config['output'] + "qc/reports/trim.p.{sample}_R1_fastqc.zip",
        config['output'] + "qc/reports/trim.p.{sample}_R2_fastqc.zip"
    threads: config["threads"]
    params:
        temp = config['output']
    shell:
        "mkdir -p {params.temp}/qc/reports; "
        "fastqc -q -t {threads} -o {params.temp}/qc/reports/ {input} "

rule map_reads:
    input:
        reference = config["reference"],
        R1 = config['output'] + "qc/data/trim.p.{sample}_R1.fastq.gz",
        R2 = config['output'] + "qc/data/trim.p.{sample}_R2.fastq.gz",
        qc_report_R1 = config['output'] + "qc/reports/trim.p.{sample}_R1_fastqc.html",
        qc_report_R2 = config['output'] + "qc/reports/trim.p.{sample}_R2_fastqc.html"
    output:
        config['output'] + "assembly/mapped_reads/{sample}.sorted.bam"
    log:
        config['output'] + "logs/minimap2/{sample}.log"
    benchmark:
        config['output'] + "logs/minimap2/{sample}.benchmark.txt"
    threads: config["threads"] 
    shell:
        "minimap2 -a -t {threads} -x sr {input.reference} {input.R1} {input.R2} | "
        "samtools view -bS -F 4 - | "
        "samtools sort -o {output} - 2> {log}"

rule index_bam_files:
    input:
        config['output'] + "assembly/mapped_reads/{sample}.sorted.bam"
    output:
        config['output'] + "assembly/mapped_reads/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule infer_consensus_sequence:
    input:
        mapped_reads = config['output'] + "assembly/mapped_reads/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/mapped_reads/{sample}.sorted.bam.bai"
    output:
        temp(config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.fasta")
    params:
        minimum_depth = config["minimum_depth"]
    benchmark:
        config['output'] + "logs/samtools/consensus/{sample}.benchmark.txt"
    shell:
        "samtools consensus -a -d {params.minimum_depth} -m simple -q -c 0.75 --show-ins yes {input.mapped_reads} -o {output}"

rule calculate_coverage_basewise:
    input:
        config['output'] + "assembly/mapped_reads/{sample}.sorted.bam"
    output:
        config['output'] + "assembly/coverage_stats/{sample}.table_cov_basewise.txt"
    shell:
        "bedtools genomecov -d -ibam {input} > {output}"

rule rename_sequences:
    input:
        config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.fasta"
    output:
        config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.renamed.fasta"
    script:
        "rename_sequences.py"

rule calculate_assembly_statistics:
    input:
        get_map_input_fastqs,
        config['output'] + "qc/data/trim.p.{sample}_R1.fastq.gz",
        config['output'] + "assembly/mapped_reads/{sample}.sorted.bam",
        config['output'] + "assembly/coverage_stats/{sample}.table_cov_basewise.txt",
        config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.renamed.fasta"
    output:
        temp(config['output'] + "assembly/coverage_stats/{sample}.stats_summary.csv")
    params:
        minimum_depth = config["minimum_depth"]
    script:
        "calculate_assembly_stats.py"

rule unify_assembly_statistics_reports:
    input:
        reports = expand(config['output'] + "assembly/coverage_stats/{sample}.stats_summary.csv", sample=config["samples"])
    output:
        config['output'] + "assembly/assembly_stats_summary.csv"
    shell:
        " echo \"sample_name,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage\" > {output} ;"
        " cat {input} >> {output}"

rule generate_multiqc_report:
    input:
        config['output'] + "assembly/assembly_stats_summary.csv"
    output:
        config['output'] + "qc/reports/multiqc_report.html"
    params:
        temp = config['output']
    shell:
        "multiqc -f -s -o {params.temp}/qc/reports/ {params.temp}/qc/reports/"

rule align_consensus_to_reference_genome:
    input:
        config['output'] + "qc/reports/multiqc_report.html"
    output:
        config['output'] + "assembly/consensus/final_consensus/aln.consensus.fasta"
    params:
        path_consensus = config['output'] + "assembly/consensus/final_consensus/",
        reference = config['reference']
    shell:
        "cat {params.reference} {params.path_consensus}/*.fasta > {params.path_consensus}/consensus.fasta; "
        "minimap2 -a --sam-hit-only --secondary=no --score-N=0 {params.reference} {params.path_consensus}/consensus.fasta -o {params.path_consensus}/aln.consensus.sam; "
        "gofasta sam toMultiAlign --pad -s {params.path_consensus}/aln.consensus.sam -o {output}"
