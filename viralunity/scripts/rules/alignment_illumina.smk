# Illumina alignment rules: perform_qc, map_reads, trim_primer_sequences
# These rules expect the following variables to be defined in the entry-point workflow:
# - REFERENCE: path to reference genome (str)
# - get_map_input_fastqs: function to resolve sample fastqs
# - config: standard Snakemake config dict


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
        reference = REFERENCE,
        R1 = rules.perform_qc.output.paired_R1,
        R2 = rules.perform_qc.output.paired_R2,
    output:
        bam = config['output'] + f"{SEGMENT_WILDCARD}assembly/mapped_reads/raw/{{sample}}.sorted.bam",
        bam_index = config['output'] + f"{SEGMENT_WILDCARD}assembly/mapped_reads/raw/{{sample}}.sorted.bam.bai"
    log:
        config['output'] + f"{SEGMENT_WILDCARD}logs/minimap2/{{sample}}.log"
    benchmark:
        config['output'] + f"{SEGMENT_WILDCARD}logs/minimap2/{{sample}}.benchmark.txt"
    threads: config["threads"] 
    shell:
        """
        minimap2 -a -t {threads} -x sr {input.reference} {input.R1} {input.R2} | \
            samtools view -bS -F 4 - | \
            samtools sort -o {output.bam} - 2> {log}
        samtools index {output.bam} {output.bam_index}
        """

rule trim_primer_sequences:
    input:
        bam = rules.map_reads.output.bam,
        bam_index = rules.map_reads.output.bam_index
    output:
        bam = config['output'] + f"{SEGMENT_WILDCARD}assembly/mapped_reads/trimmed/{{sample}}.sorted.bam",
        bam_index = config['output'] + f"{SEGMENT_WILDCARD}assembly/mapped_reads/trimmed/{{sample}}.sorted.bam.bai",
        trimmed_info = config['output'] + f"{SEGMENT_WILDCARD}assembly/mapped_reads/trimmed/{{sample}}.trimmed.txt"
    params:
        bed = config["scheme"],
        minimum_length = config["minimum_length"],
        path = config['output'] + f"{SEGMENT_WILDCARD}assembly/mapped_reads/raw/"
    log:
        config['output'] + f"{SEGMENT_WILDCARD}logs/samtools/ampliconclip/{{sample}}.log"
    benchmark:
        config['output'] + f"{SEGMENT_WILDCARD}logs/samtools/ampliconclip/{{sample}}.benchmark.txt"
    threads: config["threads"]
    shell:
        """
        if [ {params.bed} == NA ]; then
            cp {input.bam} {output.bam};
            cp {input.bam_index} {output.bam_index};
            touch {output.trimmed_info};
            echo "No primer scheme detected, assuming sequence data came from untargeted sequencing approach (bam files copied to trimmed dir for convenience)." > {params.path}notes.txt
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
