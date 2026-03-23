# Nanopore alignment rules: map_reads, trim_primer_sequences
# These rules expect the following variables to be defined in the entry-point workflow:
# - REFERENCE: path to reference genome (str)
# - config: standard Snakemake config dict


rule map_reads:
    input:
        reference = REFERENCE,
        fastq = get_map_input_fastqs
    params:
        minimum_map_quality = config["minimum_map_quality"]
    output:
        bam = config['output'] + f"{SEGMENT_WILDCARD}assembly/mapped_reads/raw/{{sample}}.sorted.bam",
        bam_index = config['output'] + f"{SEGMENT_WILDCARD}assembly/mapped_reads/raw/{{sample}}.sorted.bam.bai",
    log:
        config['output'] + f"{SEGMENT_WILDCARD}logs/minimap2/{{sample}}.log"
    benchmark:
        config['output'] + f"{SEGMENT_WILDCARD}logs/minimap2/{{sample}}.benchmark.txt"
    threads: config["threads"] 
    shell:
        """
        minimap2 -a -t {threads} -x map-ont {input.reference} {input.fastq} |
        samtools view --min-MQ {params.minimum_map_quality} -bS -F 4 - |
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
