# Shared statistics rules used by both nanopore and illumina workflows
# These rules expect the following to be defined in the entry-point workflow:
# - config: standard Snakemake config dict
# - get_map_input_fastqs: function to resolve sample fastqs
# - Rules: trim_primer_sequences, infer_consensus_sequence (from alignment/consensus includes)


rule calculate_coverage_basewise:
    conda:
        "../envs/alignment.yaml"
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        table_cov = config['output'] + "assembly/" + SEGMENT_WILDCARD + "coverage_stats/{sample}.table_cov_basewise.txt"
    shell:
        "bedtools genomecov -d -ibam {input.bam} > {output.table_cov}"


rule rename_sequences:
    conda:
        "../envs/utils.yaml"
    input:
        consensus = rules.infer_consensus_sequence.output.consensus
    output:
        consensus_renamed = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.consensus.renamed.fasta"
    script:
        "../python/rename_sequences.py"
