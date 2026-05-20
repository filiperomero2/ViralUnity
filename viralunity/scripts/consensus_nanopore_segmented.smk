SEGMENTS = config["reference"]  # dict: {"S": "/path/S.fa", "L": "/path/L.fa", ...}

rule all:
    input:
        expand(
            config['output'] + "assembly/{segment}/consensus/final_consensus/samples_alignment.fasta",
            segment=SEGMENTS.keys()
        ),
        config['output'] + "benchmark.tsv"

rule sanitize_reference:
    conda:
        "envs/clair3.yaml"
    input: lambda wildcards: SEGMENTS[wildcards.segment]
    output:
        fasta = config["output"] + "reference/{segment}.sanitized.fasta",
        fai = config["output"] + "reference/{segment}.sanitized.fasta.fai"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.fasta})
        sed '/^>/s/[\\/|,~ ]/_/g' {input} > {output.fasta}
        samtools faidx {output.fasta}
        """

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample].split()
    return(reads)

REFERENCE = rules.sanitize_reference.output.fasta
SEGMENT_WILDCARD = "{segment}/"

include: "rules/alignment_nanopore.smk"
include: "rules/consensus_nanopore.smk"
include: "rules/stats.smk"

rule calculate_assembly_statistics:
    conda:
        "envs/utils.yaml"
    input:
        get_map_input_fastqs,
        get_map_input_fastqs,
        get_map_input_fastqs,
        rules.trim_primer_sequences.output.bam,
        rules.calculate_coverage_basewise.output.table_cov,
        rules.rename_sequences.output.consensus_renamed
    output:
        stats_summary = temp(config['output'] + "assembly/{segment}/coverage_stats/{sample}.stats_summary.csv")
    params:
        minimum_depth = config["minimum_depth"]
    script:
        "python/calculate_assembly_stats.py"

rule unify_assembly_statistics_reports:
    conda:
        "envs/utils.yaml"
    input:
        reports = expand(
            rules.calculate_assembly_statistics.output.stats_summary,
            sample=config["samples"],
            segment=SEGMENTS.keys()
        )
    output:
        unified_stats_summary = config['output'] + "assembly/assembly_stats_summary.csv"
    shell:
        """
        set -euo pipefail
        echo \"sample_name,segment,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage\" > {output.unified_stats_summary} ;
        cat {input.reports} >> {output.unified_stats_summary}
        """

rule align_consensus_to_reference_genome:
    conda:
        "envs/alignment.yaml"
    input:
        stats = rules.unify_assembly_statistics_reports.output.unified_stats_summary,
        consensus_files = expand(
            rules.rename_sequences.output.consensus_renamed,
            sample=config["samples"],
            allow_missing=True
        )
    output:
        aln_consensus = config['output'] + "assembly/{segment}/consensus/final_consensus/samples_alignment.fasta"
    params:
        path_consensus = config['output'] + "assembly/{segment}/consensus/final_consensus/",
        reference = REFERENCE
    shell:
        """
        set -euo pipefail
        cat {params.reference} {params.path_consensus}/*.renamed.fasta > {params.path_consensus}/consensus.fasta;
        minimap2 -a --sam-hit-only --secondary=no --score-N=0 {params.reference} {params.path_consensus}/consensus.fasta -o {params.path_consensus}/aln.consensus.sam;
        gofasta sam toMultiAlign --pad -s {params.path_consensus}/aln.consensus.sam -o {output.aln_consensus};
        sed '/^>/ ! s/-/N/g' {output.aln_consensus} > {params.path_consensus}/aln.consensus.indelsMasked.fasta
        """

rule organize_files:
    conda:
        "envs/utils.yaml"
    input:
        vcf_files = expand(
            rules.infer_consensus_sequence.output.vcf,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        vcf_raw_files = expand(
            rules.infer_consensus_sequence.output.vcf_raw,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        table_cov = expand(
            rules.calculate_coverage_basewise.output.table_cov,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        consensus_files = expand(
            rules.rename_sequences.output.consensus_renamed,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        raw_mapped_reads = expand(
            rules.map_reads.output.bam,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        trimmed_mapped_reads = expand(
            rules.trim_primer_sequences.output.bam,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
    output:
        config['output'] + "benchmark.tsv"
    params:
        outdir = config['output'],
        samples = " ".join(config["samples"].keys()),
        segments = " ".join(SEGMENTS.keys())
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}samples/
        for sample in {params.samples}; do
            for segment in {params.segments}; do
                mkdir -p {params.outdir}samples/$sample/$segment;
            done
        done
        for _file in {input.vcf_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .vcf.gz);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/$segment/consensus.vcf.gz;
            ln -sf $PWD/$_file.tbi {params.outdir}samples/$sample/$segment/consensus.vcf.gz.tbi;
        done
        for _file in {input.vcf_raw_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo "$rel" | cut -d'/' -f1);
            sample=$(basename $_file .raw.vcf.gz);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/$segment/raw.vcf.gz;
            ln -sf $PWD/$_file.tbi {params.outdir}samples/$sample/$segment/raw.vcf.gz.tbi;
        done
        for _file in {input.table_cov}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .table_cov_basewise.txt);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/$segment/table_cov_basewise.txt;
        done
        for _file in {input.consensus_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .consensus.renamed.fasta);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/$segment/consensus.fasta;
        done
        for _file in {input.raw_mapped_reads}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .sorted.bam);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/$segment/raw_mapped_reads.bam;
            ln -sf $PWD/$_file.bai {params.outdir}samples/$sample/$segment/raw_mapped_reads.bam.bai;
        done
        for _file in {input.trimmed_mapped_reads}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .sorted.bam);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/$segment/trimmed_mapped_reads.bam;
            ln -sf $PWD/$_file.bai {params.outdir}samples/$sample/$segment/trimmed_mapped_reads.bam.bai;
        done

        # Benchmark aggregation
        echo -e "sample\\tsegment\\ttask\\tseconds\\th:m:s\\tmax_rss\\tmax_vms\\tmax_uss\\tmax_pss\\tio_in\\tio_out\\tmean_load\\tcpu_time" > {output}
        find {params.outdir} -name "*.benchmark.txt" | while read -r file; do
            task=$(basename $(dirname $file))
            sample=$(basename $file .benchmark.txt)

            outdir="{params.outdir}"; rel=${{file#$outdir}};
            if [[ "$rel" == assembly/* ]]; then
                rel=${{rel#assembly/}}
                segment=$(echo "$rel" | cut -d'/' -f1);
            else
                segment="-"
            fi

            matched=false
            for s in {params.samples}; do
                if [[ "$sample" == "$s" ]]; then
                    matched=true
                    break
                fi
            done

            if [[ "$matched" == "false" ]]; then
                sample="All"
            else
                sample=$(echo $sample | sed 's/sample-//')
            fi

            tail -n +2 $file | awk -v sample=$sample -v segment=$segment -v task=$task '{{print sample"\\t"segment"\\t"task"\\t"$0}}' >> {output}
        done
        """
