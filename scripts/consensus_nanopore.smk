# Experimental

rule all:
    input:
        config['output'] + "assembly/consensus/final_consensus/aln.consensus.fasta"

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample].split()
    return(reads)

rule map_reads:
    input:
        reference = config["reference"],
        fastq = get_map_input_fastqs
    output:
        config['output'] + "assembly/mapped_reads/raw/{sample}.sorted.bam"
    log:
        config['output'] + "logs/minimap2/{sample}.log"
    benchmark:
        config['output'] + "logs/minimap2/{sample}.benchmark.txt"
    threads: config["threads"] 
    shell:
        "minimap2 -a -t {threads} -x map-ont {input.reference} {input.fastq} | "
        "samtools view -bS -F 4 - | "
        "samtools sort -o {output} - 2> {log}"

# experimental 
rule trim_primer_sequences:
    input:
        config['output'] + "assembly/mapped_reads/raw/{sample}.sorted.bam"
    output:
        config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam"
    params:
        bed = config["scheme"],
        path = config['output'] + "assembly/mapped_reads/raw/"
    shell:
        """
        if [ {params.bed} == NA ]; then
           mv {input} {output};
           echo "No primer scheme detected, assuming sequence data came from untargeted sequencing approach (fastq files moved to trimmed dir for convenience)." > {params.path}notes.txt
        else
           ivar trim -b {params.bed} -e -q 0 -m 50 -p {params.path}trim.{wildcards.sample} -i {input};
           samtools sort {params.path}trim.{wildcards.sample}.bam -o {output};
           rm {params.path}trim.{wildcards.sample}.bam
        fi
        """

rule index_bam_files:
    input:
        config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam"
    output:
        config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule infer_consensus_sequence:
    input:
        mapped_reads = config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam.bai"
    output:
        temp(config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.fasta")
    params:
        minimum_depth = config["minimum_depth"]
    benchmark:
        config['output'] + "logs/samtools/consensus/{sample}.benchmark.txt"
    shell:
        "samtools consensus -a -d {params.minimum_depth} -m simple -q -c 0.75  -X r10.4_sup --show-ins yes {input.mapped_reads} -o {output}"


rule calculate_coverage_basewise:
    input:
        config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam"
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
        get_map_input_fastqs,
        get_map_input_fastqs,
        config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam",
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

rule align_consensus_to_reference_genome:
    input:
        config['output'] + "assembly/assembly_stats_summary.csv"
    output:
        config['output'] + "assembly/consensus/final_consensus/aln.consensus.fasta"
    params:
        path_consensus = config['output'] + "assembly/consensus/final_consensus/",
        reference = config['reference']
    shell:
        "cat {params.reference} {params.path_consensus}/*.fasta > {params.path_consensus}/consensus.fasta; "
        "minimap2 -a --sam-hit-only --secondary=no --score-N=0 {params.reference} {params.path_consensus}/consensus.fasta -o {params.path_consensus}/aln.consensus.sam; "
        "gofasta sam toMultiAlign --pad -s {params.path_consensus}/aln.consensus.sam -o {output}; "
        "sed '/^>/! s/-/N/g' {output} > {params.path_consensus}/aln.consensus.indelsMasked.fasta"
