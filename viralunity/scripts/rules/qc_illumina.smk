def get_fastp_adapter_args(wildcards):
    """Use --detect_adapter_for_pe when no adapters file is provided."""
    adapters = config.get("adapters", "NA")
    if not adapters or str(adapters).strip() in ("", "NA"):
        return "--detect_adapter_for_pe"
    return "--adapter_fasta " + str(adapters).strip()

rule perform_qc:
    conda:
        "../envs/qc.yaml"
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
        adapter_args = get_fastp_adapter_args
    threads: config.get("perform_qc_cpus", 2)
    resources:
        mem_mb = config.get("perform_qc_ram", 4) * 1024
    log:
        config['output'] + "logs/fastp/{sample}.log"
    benchmark:
        config['output'] + "logs/fastp/{sample}.benchmark.txt"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        set -euo pipefail
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
