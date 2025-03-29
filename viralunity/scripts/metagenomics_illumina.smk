rule all:
    input:
        config['output'] + "qc/reports/multiqc_report.html"

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample].split()
    return(reads)

def get_classification_filter():
    if config['remove_human_reads'] and config['remove_unclassified_reads']:
        my_filter = "\t9606|\t0"
    if not config['remove_human_reads'] and config['remove_unclassified_reads']:
        my_filter = "\t0"
    if config['remove_human_reads'] and not config['remove_unclassified_reads']:
        my_filter = "\t9606"
    if not config['remove_human_reads'] and not config['remove_unclassified_reads']:
        my_filter = " "
    print(my_filter)
    return(my_filter)

my_filter = get_classification_filter()     

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

rule run_kraken2:
    input:
        R1 = config['output'] + "qc/data/trim.p.{sample}_R1.fastq.gz",
        R2 = config['output'] + "qc/data/trim.p.{sample}_R2.fastq.gz",
        qc_report_R1 = config['output'] + "qc/reports/trim.p.{sample}_R1_fastqc.html",
        qc_report_R2 = config['output'] + "qc/reports/trim.p.{sample}_R2_fastqc.html"
    output:
        report = config['output'] + "metagenomics/taxonomic_assignments/results/{sample}.report.txt",
        outfile = config['output'] + "metagenomics/taxonomic_assignments/results/{sample}.output.txt",
    threads: 2
    params:
        database = config['kraken2_database']
    log:
        config['output'] + "logs/kraken2/{sample}.log"
    benchmark:
        config['output'] + "logs/kraken2/{sample}.benchmark.txt"
    shell:
        "kraken2 --db {params.database} --threads {threads} --report-minimizer-data "
        "--minimum-hit-group 3 --report {output.report} "
        "--output {output.outfile} --paired {input.R1} {input.R2} 2> {log}"

rule clean_kraken2_output:
    input:
        config['output'] + "metagenomics/taxonomic_assignments/results/{sample}.output.txt"
    output:
        config['output'] + "metagenomics/taxonomic_assignments/results/{sample}.output.krona.txt"
    params:
        to_rm = my_filter
    shell:
        "cat {input} | cut -f 2,3 | egrep -v \"{params.to_rm}\" > {output} "

rule create_krona_report:
    input:
        config['output'] + "metagenomics/taxonomic_assignments/results/{sample}.output.krona.txt"
    output:
        config['output'] + "metagenomics/taxonomic_assignments/reports/{sample}.output.krona.html"
    params:
        krona_database = config['krona_database']
    log:
        config['output'] + "logs/krona/{sample}.log"
    benchmark:
        config['output'] + "logs/krona/{sample}.benchmark.txt"
    shell:
        "ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}"

rule summarize:
    input:
        expand(config['output'] + "metagenomics/taxonomic_assignments/reports/{sample}.output.krona.html", sample=config["samples"])
    output:
        config['output'] + "metagenomics/metagenomics_summary.txt"
    params:
        path = config['output'] + "metagenomics/taxonomic_assignments/results/"
    shell:
        "echo \"file,percentage_of_reads,number_of_reads_rooted,number_of_reads_direct,number_of_k-mers,number_of_distinct_k-mers,rank_code,ncbi_taxid,taxon\" > {output}; "
        "grep \"viridae\" {params.path}*report.txt | sed -E \"s/[	| |:]+/,/g\">> {output}"

rule generate_multiqc_report:
    input:
        config['output'] + "metagenomics/metagenomics_summary.txt"
    output:
        config['output'] + "qc/reports/multiqc_report.html"
    params:
        temp = config['output']
    shell:
        "multiqc -f -s -o {params.temp}/qc/reports/ {params.temp}/qc/reports/"
