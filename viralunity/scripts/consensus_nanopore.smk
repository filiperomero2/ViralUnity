SEGMENT_WILDCARD = ""
REFERENCE = config["reference"]

rule all:
    input:
        config['output'] + "assembly/consensus/final_consensus/aln.consensus.fasta",
        config['output'] + "benchmark.tsv"

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample].split()
    return(reads)

include: "rules/alignment_nanopore.smk"
include: "rules/consensus_nanopore.smk"
include: "rules/stats.smk"


# The calculate_assembly_stats.py expects three fastq inputs: raw_r1, raw_r2 and trimmed
# Here we employ a workaround by passing the same fastq file for raw_r1 and raw_r2
# in Nanopore scenario
rule calculate_assembly_statistics:
    input:
        get_map_input_fastqs,
        get_map_input_fastqs,
        get_map_input_fastqs,
        rules.trim_primer_sequences.output.bam,
        rules.calculate_coverage_basewise.output.table_cov,
        rules.rename_sequences.output.consensus_renamed
    output:
        stats_summary = temp(config['output'] + "assembly/coverage_stats/{sample}.stats_summary.csv")
    params:
        minimum_depth = config["minimum_depth"]
    script:
        "calculate_assembly_stats.py"

rule unify_assembly_statistics_reports:
    input:
        reports = expand(rules.calculate_assembly_statistics.output.stats_summary, sample=config["samples"])
    output:
        unified_stats_summary = config['output'] + "assembly/assembly_stats_summary.csv"
    shell:
        """
        echo \"sample_name,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage\" > {output.unified_stats_summary} ;
        cat {input.reports} >> {output.unified_stats_summary}
        """


rule align_consensus_to_reference_genome:
    input:
        rules.unify_assembly_statistics_reports.output.unified_stats_summary
    output:
        aln_consensus = config['output'] + "assembly/consensus/final_consensus/aln.consensus.fasta"
    params:
        path_consensus = config['output'] + "assembly/consensus/final_consensus/",
        reference = REFERENCE
    shell:
        """
        cat {params.reference} {params.path_consensus}/*.renamed.fasta > {params.path_consensus}/consensus.fasta; 
        minimap2 -a --sam-hit-only --secondary=no --score-N=0 {params.reference} {params.path_consensus}/consensus.fasta -o {params.path_consensus}/aln.consensus.sam; 
        gofasta sam toMultiAlign --pad -s {params.path_consensus}/aln.consensus.sam -o {output.aln_consensus}; 
        sed '/^>/ ! s/-/N/g' {output.aln_consensus} > {params.path_consensus}/aln.consensus.indelsMasked.fasta
        """

rule organize_files:
    input:
        vcf_files = expand(rules.infer_consensus_sequence.output.vcf, sample=config["samples"]),
        vcf_raw_files = expand(rules.infer_consensus_sequence.output.vcf_raw, sample=config["samples"]),
        table_cov = expand(rules.calculate_coverage_basewise.output.table_cov, sample=config["samples"]),
        consensus_files = expand(rules.rename_sequences.output.consensus_renamed, sample=config["samples"]),
        raw_mapped_reads = expand(rules.map_reads.output.bam, sample=config["samples"]),
        trimmed_mapped_reads = expand(rules.trim_primer_sequences.output.bam, sample=config["samples"]),
    output:
        config['output'] + "benchmark.tsv"
    params:
        outdir = config['output'],
        samples = " ".join(config["samples"].keys())
    shell:
        """
        mkdir -p {params.outdir}samples/
        for sample in {params.samples}; do
            mkdir -p {params.outdir}samples/$sample;
        done
        for _file in {input.vcf_files}; do
            sample=$(basename $_file .vcf.gz);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/consensus.vcf.gz;
            ln -sf $PWD/$_file.tbi {params.outdir}samples/$sample/consensus.vcf.gz.tbi;
        done
        for _file in {input.vcf_raw_files}; do
            sample=$(basename $_file .raw.vcf.gz);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/raw.vcf.gz;
            ln -sf $PWD/$_file.tbi {params.outdir}samples/$sample/raw.vcf.gz.tbi;
        done
        for _file in {input.table_cov}; do
            sample=$(basename $_file .table_cov_basewise.txt);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/table_cov_basewise.txt;
        done
        for _file in {input.consensus_files}; do
            sample=$(basename $_file .consensus.renamed.fasta);
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