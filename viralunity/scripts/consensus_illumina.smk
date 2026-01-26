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
        paired_R2 = config['output'] + "qc/data/trim.p.{sample}_R2.fastq.gz",
        unpaired_R1 = config['output'] + "qc/data/trim.u.{sample}_R1.fastq.gz",
        unpaired_R2 = config['output'] + "qc/data/trim.u.{sample}_R2.fastq.gz",
        json = config['output'] + "qc/reports/trim.{sample}_fastp.json",
        html = config['output'] + "qc/reports/trim.{sample}_fastp.html"
    params:
        minimum_length = config["minimum_length"],
        adapters = config["adapters"],
        trim_head = config["trim_head"],
        trim_tail = config["trim_tail"],
        cut_front_mean_quality = config["cut_front_mean_quality"],
        cut_tail_mean_quality = config["cut_tail_mean_quality"],
        cut_right_window_size = config["cut_right_window_size"],
        cut_right_mean_quality = config["cut_right_mean_quality"],
        adapter_args = lambda w: "--detect_adapter_for_pe" if not config["adapters"] else "--adapter_fasta " + config["adapters"]
    threads: config["threads"]
    log:
        config['output'] + "logs/fastp/{sample}.log"
    benchmark:
        config['output'] + "logs/fastp/{sample}.benchmark.txt"
    shell:
        """
        fastp \
            -i {input[0]} \
            -I {input[1]} \
            -o {output.paired_R1} \
            -O {output.paired_R2} \
            --unpaired1 {output.unpaired_R1} \
            --unpaired2 {output.unpaired_R2} \
            --length_required {params.minimum_length} \
            --trim_front1 {params.trim_head} \
            --trim_tail1 {params.trim_tail} \
            --trim_front2 {params.trim_head} \
            --trim_tail2 {params.trim_tail} \
            --cut_front \
            --cut_front_mean_quality {params.cut_front_mean_quality} \
            --cut_tail \
            --cut_tail_mean_quality {params.cut_tail_mean_quality} \
            --cut_right \
            --cut_right_window_size {params.cut_right_window_size} \
            --cut_right_mean_quality {params.cut_right_mean_quality} \
            --json {output.json} \
            --html {output.html} \
            --thread {threads} \
            {params.adapter_args}
        """

rule map_reads:
    input:
        reference = config["reference"],
        R1 = rules.perform_qc.output.paired_R1,
        R2 = rules.perform_qc.output.paired_R2,
    output:
        bam = config['output'] + "assembly/mapped_reads/raw/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/mapped_reads/raw/{sample}.sorted.bam.bai"
    log:
        config['output'] + "logs/minimap2/{sample}.log"
    benchmark:
        config['output'] + "logs/minimap2/{sample}.benchmark.txt"
    threads: config["threads"] 
    shell:
        """
        minimap2 -a -t {threads} -x sr {input.reference} {input.R1} {input.R2} | \
            samtools view -bS -F 4 - | \
            samtools sort -o {output.bam} - 2> {log}
        samtools index {output.bam} {output.bam_index}
        """

# experimental 
rule trim_primer_sequences:
    input:
        bam = rules.map_reads.output.bam,
        bam_index = rules.map_reads.output.bam_index
    output:
        bam = config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam.bai"
    params:
        bed = config["scheme"],
        path = config['output'] + "assembly/mapped_reads/raw/"
    shell:
        """
        if [ {params.bed} == NA ]; then
           mv {input.bam} {output.bam};
           mv {input.bam_index} {output.bam_index};
           echo "No primer scheme detected, assuming sequence data came from untargeted sequencing approach (fastq files moved to trimmed dir for convenience)." > {params.path}notes.txt
        else
           ivar trim -b {params.bed} -e -q 0 -m 50 -p {params.path}trim.{wildcards.sample} -i {input.bam};
           samtools sort {params.path}trim.{wildcards.sample}.bam -o {output.bam};
           rm {params.path}trim.{wildcards.sample}.bam
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

rule calculate_assembly_statistics:
    input:
        get_map_input_fastqs,
        rules.perform_qc.output.paired_R1,
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

rule generate_multiqc_report:
    input:
        unified_stats_summary = rules.unify_assembly_statistics_reports.output.unified_stats_summary
    output:
        multiqc_report = config['output'] + "qc/reports/multiqc_report.html"
    params:
        temp = config['output']
    shell:
        "multiqc -f -o {params.temp}/qc/reports/ {params.temp}/qc/reports/"

rule align_consensus_to_reference_genome:
    input:
        rules.generate_multiqc_report.output.multiqc_report
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
