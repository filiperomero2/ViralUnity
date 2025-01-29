rule all:
    input:
        config['output'] + "qc/reports/multiqc_report.html"

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

rule generate_multiqc_report:
    input:
        expand(config['output'] + "qc/reports/trim.p.{sample}_R1_fastqc.html", sample=config["samples"]),
        expand(config['output'] + "qc/reports/trim.p.{sample}_R2_fastqc.html", sample=config["samples"])
    output:
        config['output'] + "qc/reports/multiqc_report.html"
    params:
        temp = config['output']
    shell:
        "multiqc -f -s -o {params.temp}/qc/reports/ {params.temp}/qc/reports/"
