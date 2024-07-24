rule all:
    input:
        config['output'] + "assembly/annotation/mutation_report.txt"

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

rule call_variants:
    input:
        reference = config["reference"],
        mapped_reads = config['output'] + "assembly/mapped_reads/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/mapped_reads/{sample}.sorted.bam.bai"
    output:
        config['output'] + "assembly/variant_calls/{sample}.calls.vcf.gz"
    threads: config["threads"]
    log:
        config['output'] + "logs/bcftools/raw_calls/{sample}.log"
    benchmark:
        config['output'] + "logs/bcftools/raw_calls/{sample}.benchmark.txt"
    shell:
        "bcftools mpileup --threads {threads} "
        "--max-depth 200 -E -Ou "
        "-f {input.reference} {input.mapped_reads}  | "
        "bcftools call --ploidy 1 --threads {threads} "
        " -mv -Oz -o {output} 2> {log}"

rule normalize_variants:
    input:
        reference = config["reference"],
        variants = config['output'] + "assembly/variant_calls/{sample}.calls.vcf.gz"
    output:
        config['output'] + "assembly/variant_calls/{sample}.calls.norm.vcf.gz"
    threads: config["threads"]
    log:
        config['output'] + "logs/bcftools/normalized_calls/{sample}.log"
    benchmark:
        config['output'] + "logs/bcftools/normalized_calls/{sample}.benchmark.txt"
    
    shell:
        "bcftools norm --threads {threads} -f {input.reference} "
        "{input.variants} -Oz -o {output} 2> {log}"

rule filter_normalized_variants:
    input:
        config['output'] + "assembly/variant_calls/{sample}.calls.norm.vcf.gz"
    output:
        config['output'] + "assembly/variant_calls/{sample}.calls.norm.flt-indels.vcf.gz"
    threads: config["threads"]
    log:
        config['output'] + "logs/bcftools/filtered_calls/{sample}.log"
    benchmark:
        config['output'] + "logs/bcftools/filtered_calls/{sample}.benchmark.txt"
    shell:
        "bcftools filter --threads {threads} "
        "--IndelGap 5 {input} "
        "-Oz -o {output}"

rule index_normalized_filtered_variants:
    input:
        config['output'] + "assembly/variant_calls/{sample}.calls.norm.flt-indels.vcf.gz"
    output:
        config['output'] + "assembly/variant_calls/{sample}.calls.norm.flt-indels.vcf.gz.csi"
    threads: config["threads"]
    log:
        config['output'] + "logs/bcftools/indexed_calls/{sample}.log"
    benchmark:
        config['output'] + "logs/bcftools/indexed_calls/{sample}.benchmark.txt"
    shell:
        "bcftools index --threads {threads} {input} 2> {log}"

rule infer_raw_consensus_sequence:
    input:
        reference = config["reference"],
        variants = config['output'] + "assembly/variant_calls/{sample}.calls.norm.flt-indels.vcf.gz",
        variants_index = config['output'] + "assembly/variant_calls/{sample}.calls.norm.flt-indels.vcf.gz.csi"
    output:
        temp(config['output'] + "assembly/consensus/raw_consensus/{sample}.temp.consensus.fasta")
    log:
        config['output'] + "logs/bcftools/consensus/{sample}.log"
    benchmark:
        config['output'] + "logs/bcftools/consensus/{sample}.benchmark.txt"
    shell:
        "bcftools consensus -f {input.reference} "
        "{input.variants} > {output} 2> {log} "

rule calculate_coverage:
    input:
        config['output'] + "assembly/mapped_reads/{sample}.sorted.bam"
    output:
        config['output'] + "assembly/coverage_stats/{sample}.table_cov.txt"
    params:
        minimum_depth = config["minimum_depth"]
    shell:
        "bedtools genomecov -bga -ibam {input} > {output}"

rule calculate_coverage_basewise:
    input:
        config['output'] + "assembly/mapped_reads/{sample}.sorted.bam"
    output:
        config['output'] + "assembly/coverage_stats/{sample}.table_cov_basewise.txt"
    shell:
        "bedtools genomecov -d -ibam {input} > {output}"

rule create_minimum_depth_table:
    input:
        config['output'] + "assembly/coverage_stats/{sample}.table_cov.txt"
    output:
        config['output'] + "assembly/coverage_stats/{sample}.table_cov_minimum_depth.txt"
    params:
        minimum_depth = config["minimum_depth"]
    shell:
        "awk  '$4 < '{params.minimum_depth}'' {input} > {output}"

rule mask_raw_consensus_sequence:
    input:
        raw_consensus = config['output'] + "assembly/consensus/raw_consensus/{sample}.temp.consensus.fasta",
        minimum_depth_table = config['output'] + "assembly/coverage_stats/{sample}.table_cov_minimum_depth.txt"
    output:
        temp(config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.fasta")
    script:
        "mask_consensus.py"

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

rule create_mutations_report:
    input:
        consensus = config['output'] + "assembly/consensus/final_consensus/aln.consensus.fasta",
        annotation = config['annotation']
    output:
        config['output'] + "assembly/annotation/mutation_report.txt"
    params:
        path = config['workflow_path']
    shell:
        "python {params.path}/mutation_mapper.py --input {input.consensus} --annotation-file {input.annotation} --output {output}"
      