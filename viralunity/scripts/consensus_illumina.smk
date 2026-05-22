SEGMENT_WILDCARD = ""
REFERENCE = config["reference"]
SAMPLE_SEGMENT_LOG_PREFIX = config['output'] + "logs/consensus_illumina/"

rule all:
    input:
        config['output'] + "assembly/consensus/final_consensus/samples_alignment.fasta",
        config['output'] + "isnvs/isnvs_summary.tsv" if config.get("run_isnv", False) else [],
        config['output'] + "benchmark.tsv"

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample].split()
    return(reads)

include: "rules/qc_illumina.smk"
include: "rules/alignment_illumina.smk"
include: "rules/consensus_illumina.smk"
include: "rules/stats.smk"
include: "rules/consensus_illumina_common.smk"

rule unify_assembly_statistics_reports:
    conda:
        "envs/utils.yaml"
    input:
        reports = expand(rules.calculate_assembly_statistics.output.stats_summary, sample=config["samples"])
    output:
        unified_stats_summary = config['output'] + "assembly/assembly_stats_summary.csv"
    log:
        config['output'] + "logs/consensus_illumina/unify_assembly_statistics_reports/unify_assembly_statistics_reports.log"
    benchmark:
        config['output'] + "logs/consensus_illumina/unify_assembly_statistics_reports/unify_assembly_statistics_reports.benchmark.txt"
    shell:
        """
        set -euo pipefail
        echo \"sample_name,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage\" > {output.unified_stats_summary} ;
        cat {input.reports} >> {output.unified_stats_summary}
        """

rule summarize_isnvs:
    conda:
        "envs/utils.yaml"
    input:
        vcf_files = expand(rules.detect_isnv.output.vcf, sample=config["samples"])
    output:
        isnvs_summary = config['output'] +  "isnvs/isnvs_summary.tsv"
    log:
        config['output'] + "logs/consensus_illumina/summarize_isnvs/summarize_isnvs.log"
    benchmark:
        config['output'] + "logs/consensus_illumina/summarize_isnvs/summarize_isnvs.benchmark.txt"   
    shell:
        """
        set -euo pipefail
        echo -e "sample\\tnumber_of_isnvs" > {output.isnvs_summary};
        for _file in {input.vcf_files}; do
            sample=$(basename $_file .isnvs.vcf.gz);
            isnv_count=$(bcftools view -H $_file | wc -l);
            echo -e "$sample\\t$isnv_count" >> {output.isnvs_summary};
        done
        """

rule organize_files:
    conda:
        "envs/utils.yaml"
    input:
        fastp_reports = expand(rules.perform_qc.output.html, sample=config["samples"]),
        vcf_files = expand(rules.generate_vcf_consensus.output.vcf, sample=config["samples"]),
        isn_vcf_files = expand(rules.detect_isnv.output.vcf, sample=config["samples"]) if config.get("run_isnv", False) else [],
        stats_summary = expand(rules.calculate_assembly_statistics.output.stats_summary, sample=config["samples"]),
        consensus_files = expand(rules.infer_consensus_sequence.output.consensus, sample=config["samples"]),
        raw_mapped_reads = expand(rules.map_reads.output.bam, sample=config["samples"]),
        trimmed_mapped_reads = expand(rules.trim_primer_sequences.output.bam, sample=config["samples"]),
    output:
        config['output'] + "benchmark.tsv"
    params:
        outdir = config['output'],
        samples = " ".join(config["samples"].keys())
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}samples/
        for sample in {params.samples}; do
            mkdir -p {params.outdir}samples/$sample;
        done
        for _file in {input.fastp_reports}; do
            sample=$(basename $_file _fastp.html | sed 's/^trim.//');
            ln -sf $PWD/$_file {params.outdir}samples/$sample/fastp.html;
        done
        for _file in {input.vcf_files}; do
            sample=$(basename $_file .consensus.vcf.gz);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/consensus.vcf.gz;
            ln -sf $PWD/$_file.tbi {params.outdir}samples/$sample/consensus.vcf.gz.tbi;
        done
        for _file in {input.isn_vcf_files} ""; do
            if [ -z "$_file" ]; then continue; fi
            sample=$(basename $_file .isnvs.vcf.gz);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/isnvs.vcf.gz;
            ln -sf $PWD/$_file.tbi {params.outdir}samples/$sample/isnvs.vcf.gz.tbi;
        done
        for _file in {input.stats_summary}; do
            sample=$(basename $_file .stats_summary.csv);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/stats_summary.csv;
        done
        for _file in {input.consensus_files}; do
            sample=$(basename $_file .consensus.fasta);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/consensus.fasta;
        done
        for _file in {input.raw_mapped_reads}; do
            sample=$(basename $_file .sorted.bam);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/raw_mapped_reads.bam;
            ln -sf $PWD/$_file.bai {params.outdir}samples/$sample/raw_mapped_reads.bam.bai;
        done
        for _file in {input.trimmed_mapped_reads}; do
            sample=$(basename $_file .sorted.bam);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/trimmed_mapped_reads.bam;
            ln -sf $PWD/$_file.bai {params.outdir}samples/$sample/trimmed_mapped_reads.bam.bai;
        done
        
        # Benchmark aggregation
        echo -e "sample\\ttask\\tseconds\\th:m:s\\tmax_rss\\tmax_vms\\tmax_uss\\tmax_pss\\tio_in\\tio_out\\tmean_load\\tcpu_time" > {output}
        find {params.outdir} -name "*.benchmark.txt" | while read -r file; do
            task=$(basename $(dirname $file))
            sample=$(basename $file .benchmark.txt)
            
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

            tail -n +2 $file | awk -v sample=$sample -v task=$task '{{print sample"\\t"task"\\t"$0}}' >> {output}
        done
        """
