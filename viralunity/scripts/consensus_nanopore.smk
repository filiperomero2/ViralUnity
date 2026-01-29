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
        bam = config['output'] + "assembly/mapped_reads/raw/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/mapped_reads/raw/{sample}.sorted.bam.bai",
    log:
        config['output'] + "logs/minimap2/{sample}.log"
    benchmark:
        config['output'] + "logs/minimap2/{sample}.benchmark.txt"
    threads: config["threads"] 
    shell:
        """
        minimap2 -a -t {threads} -x map-ont {input.reference} {input.fastq} |
        samtools view -bS -F 4 - |
        samtools sort -o {output.bam} - 2> {log}
        samtools index {output.bam} {output.bam_index}
        """

rule trim_primer_sequences:
    input:
        bam = rules.map_reads.output.bam,
        bam_index = rules.map_reads.output.bam_index
    output:
        bam = config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam.bai",
        trimmed_info = config['output'] + "assembly/mapped_reads/trimmed/{sample}.trimmed.txt"
    params:
        bed = config["scheme"],
        minimum_length = config["minimum_length"],
        path = config['output'] + "assembly/mapped_reads/raw/"
    shell:
        """
        if [ {params.bed} == NA ]; then
            mv {input.bam} {output.bam};
            mv {input.bam_index} {output.bam_index};
            touch {output.trimmed_info};
            echo "No primer scheme detected, assuming sequence data came from untargeted sequencing approach (fastq files moved to trimmed dir for convenience)." > {params.path}notes.txt
        else
            samtools ampliconclip \
                --both-ends \
                --hard-clip \
                --filter-len {params.minimum_length} \
                -b {params.bed} \
                -f {output.trimmed_info} \
                {input.bam} > {output.bam}
            samtools index {output.bam} {output.bam_index}
        fi
        """

rule infer_consensus_sequence:
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        consensus = temp(config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.fasta")
    params:
        minimum_depth = config["minimum_depth"],
        af_threshold = config["af_threshold"]
    benchmark:
        config['output'] + "logs/samtools/consensus/{sample}.benchmark.txt"
    shell:
        "samtools consensus -a -d {params.minimum_depth} -m simple -q -c {params.af_threshold} --show-ins yes {input.bam} -o {output.consensus}"



rule calculate_coverage_basewise:
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        table_cov = config['output'] + "assembly/coverage_stats/{sample}.table_cov_basewise.txt"
    shell:
        "bedtools genomecov -d -ibam {input.bam} > {output.table_cov}"


rule rename_sequences:
    input:
        consensus = rules.infer_consensus_sequence.output.consensus
    output:
        consensus_renamed = config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.renamed.fasta"
    script:
        "rename_sequences.py"

# The calculate_assembly_stats.py expects three fastq inputs: raw_r1, raw_r2 and trimmed
# Here we employ a workaround by passing the same fastq file for raw_r1 and raw_r2
# in Nanopore scenario
rule calculate_assembly_statistics:
    input:
        get_map_input_fastqs,
        get_map_input_fastqs,
        get_map_input_fastqs,
        rules.trim_primer_sequences.output.bam,
        rules.calculate_coverage_basewise.output.table_cov,
        rules.rename_sequences.output.consensus_renamed
    output:
        stats_summary = temp(config['output'] + "assembly/coverage_stats/{sample}.stats_summary.csv")
    params:
        minimum_depth = config["minimum_depth"]
    script:
        "calculate_assembly_stats.py"

rule unify_assembly_statistics_reports:
    input:
        reports = expand(rules.calculate_assembly_statistics.output.stats_summary, sample=config["samples"])
    output:
        unified_stats_summary = config['output'] + "assembly/assembly_stats_summary.csv"
    shell:
        """
        echo \"sample_name,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage\" > {output.unified_stats_summary} ;
        cat {input.reports} >> {output.unified_stats_summary}
        """


rule align_consensus_to_reference_genome:
    input:
        rules.unify_assembly_statistics_reports.output.unified_stats_summary
    output:
        aln_consensus = config['output'] + "assembly/consensus/final_consensus/aln.consensus.fasta"
    params:
        path_consensus = config['output'] + "assembly/consensus/final_consensus/",
        reference = config['reference']
    shell:
        """
        cat {params.reference} {params.path_consensus}/*.fasta > {params.path_consensus}/consensus.fasta; 
        minimap2 -a --sam-hit-only --secondary=no --score-N=0 {params.reference} {params.path_consensus}/consensus.fasta -o {params.path_consensus}/aln.consensus.sam; 
        gofasta sam toMultiAlign --pad -s {params.path_consensus}/aln.consensus.sam -o {output.aln_consensus}; 
        sed '/^>/! s/-/N/g' {output.aln_consensus} > {params.path_consensus}/aln.consensus.indelsMasked.fasta
        """