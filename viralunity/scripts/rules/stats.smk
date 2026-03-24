# Shared statistics rules used by both nanopore and illumina workflows
# These rules expect the following to be defined in the entry-point workflow:
# - config: standard Snakemake config dict
# - get_map_input_fastqs: function to resolve sample fastqs
# - Rules: trim_primer_sequences, infer_consensus_sequence (from alignment/consensus includes)


rule calculate_coverage_basewise:
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        table_cov = config['output'] + f"assembly/{SEGMENT_WILDCARD}coverage_stats/{{sample}}.table_cov_basewise.txt"
    shell:
        "bedtools genomecov -d -ibam {input.bam} > {output.table_cov}"


rule rename_sequences:
    input:
        consensus = rules.infer_consensus_sequence.output.consensus
    output:
        consensus_renamed = config['output'] + f"assembly/{SEGMENT_WILDCARD}consensus/final_consensus/{{sample}}.consensus.renamed.fasta"
    script:
        "../rename_sequences.py"
