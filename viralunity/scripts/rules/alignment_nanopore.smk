# Nanopore alignment rules: map_reads, trim_primer_sequences
# These rules expect the following variables to be defined in the entry-point workflow:
# - REFERENCE: path to reference genome (str)
# - config: standard Snakemake config dict


rule map_reads:
    conda:
        "../envs/alignment.yaml"
    input:
        reference = REFERENCE,
        fastq = get_map_input_fastqs
    params:
        minimum_map_quality = config.get("minimum_map_quality", 20)
    output:
        bam = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/raw/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/raw/{sample}.sorted.bam.bai",
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/minimap2/{sample}.log"
    benchmark:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/minimap2/{sample}.benchmark.txt"
    threads: config.get("map_reads_cpus", 2)
    resources:
        mem_mb = config.get("map_reads_ram", 4) * 1024
    shell:
        """
        minimap2 -a -t {threads} -x map-ont {input.reference} {input.fastq} |
        samtools view --min-MQ {params.minimum_map_quality} -bS -F 4 - |
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
        bed = config.get("scheme", "NA"),
        minimum_length = config.get("minimum_length", 200),
        path = config['output'] + "assembly/" + SEGMENT_WILDCARD + "mapped_reads/raw/"
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
