# Illumina alignment rules: map_reads, trim_primer_sequences
# These rules expect the following variables to be defined in the entry-point workflow:
# - REFERENCE: path to reference genome (str)
# - get_map_input_fastqs: function to resolve sample fastqs
# - config: standard Snakemake config dict


rule map_reads:
    conda:
        "../envs/alignment.yaml"
    input:
        reference = REFERENCE,
        R1 = rules.perform_qc.output.paired_R1,
        R2 = rules.perform_qc.output.paired_R2,
    output:
        bam = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/raw/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/raw/{sample}.sorted.bam.bai"
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/minimap2/{sample}.log"
    benchmark:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/minimap2/{sample}.benchmark.txt"
    threads: config.get("map_reads_cpus", 2)
    resources:
        mem_mb = config.get("map_reads_ram", 4) * 1024
    shell:
        """
        minimap2 -a -t {threads} -x sr {input.reference} {input.R1} {input.R2} | \
            samtools view -bS -F 4 - | \
            samtools sort -o {output.bam} - 2> {log}
        samtools index {output.bam} {output.bam_index}
        """

rule trim_primer_sequences:
    conda:
        "../envs/alignment.yaml"
    input:
        bam = rules.map_reads.output.bam,
        bam_index = rules.map_reads.output.bam_index
    output:
        bam = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/trimmed/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/trimmed/{sample}.sorted.bam.bai",
        trimmed_info = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/trimmed/{sample}.trimmed.txt"
    params:
        bed = config["scheme"],
        minimum_length = config["minimum_length"],
        path = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/raw/"
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/samtools/ampliconclip/{sample}.log"
    benchmark:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/samtools/ampliconclip/{sample}.benchmark.txt"
    threads: config.get("trim_primer_sequences_cpus", 2)
    resources:
        mem_mb = config.get("trim_primer_sequences_ram", 4) * 1024
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
