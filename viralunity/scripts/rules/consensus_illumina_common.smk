"""Rules shared between the segmented and non-segmented Illumina consensus
workflows.

The including snakefile must define:
* ``SEGMENT_WILDCARD``       -- ``"{segment}/"`` (segmented) or ``""`` (single)
* ``REFERENCE``              -- path string (single) or callable taking
                                ``wildcards`` and returning a path (segmented)
* ``SAMPLE_SEGMENT_LOG_PREFIX`` -- log/benchmark prefix for per-sample (and,
                                in segmented mode, per-segment) rules
"""

_TOP_LOG = config['output'] + "logs/consensus_illumina/"


rule calculate_assembly_statistics:
    conda:
        "../envs/utils.yaml"
    input:
        get_map_input_fastqs,
        rules.perform_qc.output.paired_R1,
        rules.trim_primer_sequences.output.bam,
        rules.calculate_coverage_basewise.output.table_cov,
        rules.rename_sequences.output.consensus_renamed
    output:
        stats_summary = config['output'] + "assembly/" + SEGMENT_WILDCARD + "coverage_stats/{sample}.stats_summary.csv"
    params:
        minimum_depth = config["minimum_depth"]
    log:
        SAMPLE_SEGMENT_LOG_PREFIX + "calculate_assembly_statistics/{sample}.log"
    benchmark:
        SAMPLE_SEGMENT_LOG_PREFIX + "calculate_assembly_statistics/{sample}.benchmark.txt"
    script:
        "../python/calculate_assembly_stats.py"


rule generate_multiqc_report:
    conda:
        "../envs/qc.yaml"
    input:
        # Literal path (instead of ``rules.unify_assembly_statistics_reports.``)
        # to avoid an include-order cycle: ``unify_assembly_statistics_reports``
        # lives in the top-level snakefile and itself references
        # ``rules.calculate_assembly_statistics`` defined here.
        unified_stats_summary = config['output'] + "assembly/assembly_stats_summary.csv"
    output:
        multiqc_report = config['output'] + "qc/reports/multiqc_report.html"
    params:
        temp = config['output']
    log:
        _TOP_LOG + "generate_multiqc_report/generate_multiqc_report.log"
    benchmark:
        _TOP_LOG + "generate_multiqc_report/generate_multiqc_report.benchmark.txt"
    shell:
        "multiqc -f -o {params.temp}/qc/reports/ {params.temp}/qc/reports/"


rule align_consensus_to_reference_genome:
    conda:
        "../envs/alignment.yaml"
    input:
        rules.generate_multiqc_report.output.multiqc_report,
        consensus_files = expand(
            rules.rename_sequences.output.consensus_renamed,
            sample=config["samples"],
            allow_missing=True
        )
    output:
        aln_consensus = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/samples_alignment.fasta"
    params:
        path_consensus = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/",
        reference = REFERENCE,
        minimap2_flags = config.get("minimap2_consensus_align_flags", "-a --sam-hit-only --secondary=no --score-N=0")
    log:
        SAMPLE_SEGMENT_LOG_PREFIX + "align_consensus_to_reference_genome/align_consensus_to_reference_genome.log"
    benchmark:
        SAMPLE_SEGMENT_LOG_PREFIX + "align_consensus_to_reference_genome/align_consensus_to_reference_genome.benchmark.txt"
    shell:
        """
        set -euo pipefail
        cat {params.reference} {params.path_consensus}/*renamed.fasta > {params.path_consensus}/consensus.fasta;
        minimap2 {params.minimap2_flags} {params.reference} {params.path_consensus}/consensus.fasta -o {params.path_consensus}/aln.consensus.sam;
        gofasta sam toMultiAlign --pad -s {params.path_consensus}/aln.consensus.sam -o {output.aln_consensus};
        sed '/^>/ ! s/-/N/g' {output.aln_consensus} > {params.path_consensus}/aln.consensus.indelsMasked.fasta
        """
