SEGMENTS = config["reference"]  # dict: {"S": "/path/S.fa", "L": "/path/L.fa", ...}

rule all:
    input:
        expand(
            config['output'] + "assembly/{segment}/consensus/final_consensus/aln.consensus.fasta",
            segment=SEGMENTS.keys()
        ),
        config['output'] + "isnvs/isnvs_summary.tsv" if config.get("run_isnv", False) else [],
        config['output'] + "benchmark.tsv"

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample].split()
    return(reads)

def get_segment_reference(wildcards):
    return SEGMENTS[wildcards.segment]

REFERENCE = get_segment_reference
SEGMENT_WILDCARD = "{segment}/"

include: "rules/alignment_illumina.smk"
include: "rules/consensus_illumina.smk"
include: "rules/stats.smk"

rule calculate_assembly_statistics:
    input:
        get_map_input_fastqs,
        rules.perform_qc.output.paired_R1,
        rules.trim_primer_sequences.output.bam,
        rules.calculate_coverage_basewise.output.table_cov,
        rules.rename_sequences.output.consensus_renamed
    output:
        stats_summary = config['output'] + "assembly/{segment}/coverage_stats/{sample}.stats_summary.csv"
    params:
        minimum_depth = config["minimum_depth"]
    log:
        config['output'] + "assembly/{segment}/logs/consensus_illumina/calculate_assembly_statistics/{sample}.log"
    benchmark:
        config['output'] + "assembly/{segment}/logs/consensus_illumina/calculate_assembly_statistics/{sample}.benchmark.txt"
    script:
        "calculate_assembly_stats.py"

rule unify_assembly_statistics_reports:
    input:
        reports = expand(
            rules.calculate_assembly_statistics.output.stats_summary,
            sample=config["samples"],
            segment=SEGMENTS.keys()
        )
    output:
        unified_stats_summary = config['output'] + "assembly/assembly_stats_summary.csv"
    log:
        config['output'] + "logs/consensus_illumina/unify_assembly_statistics_reports/unify_assembly_statistics_reports.log"
    benchmark:
        config['output'] + "logs/consensus_illumina/unify_assembly_statistics_reports/unify_assembly_statistics_reports.benchmark.txt"
    shell:
        """
        echo \"sample_name,segment,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage\" > {output.unified_stats_summary} ;
        cat {input.reports} >> {output.unified_stats_summary}
        """

rule generate_multiqc_report:
    input:
        unified_stats_summary = rules.unify_assembly_statistics_reports.output.unified_stats_summary
    output:
        multiqc_report = config['output'] + "qc/reports/multiqc_report.html"
    params:
        temp = config['output']
    log:
        config['output'] + "logs/consensus_illumina/generate_multiqc_report/generate_multiqc_report.log"
    benchmark:
        config['output'] + "logs/consensus_illumina/generate_multiqc_report/generate_multiqc_report.benchmark.txt"
    shell:
        "multiqc -f -o {params.temp}/qc/reports/ {params.temp}/qc/reports/"

rule summarize_isnvs:
    input:
        vcf_files = expand(
            rules.detect_isnv.output.vcf,
            sample=config["samples"],
            segment=SEGMENTS.keys()
        ) if config.get("run_isnv", False) else []
    output:
        isnvs_summary = config['output'] + "isnvs/isnvs_summary.tsv"
    log:
        config['output'] + "logs/consensus_illumina/summarize_isnvs/summarize_isnvs.log"
    benchmark:
        config['output'] + "logs/consensus_illumina/summarize_isnvs/summarize_isnvs.benchmark.txt"
    params:
        outdir = config['output']
    shell:
        """
        echo -e "sample\\tsegment\\tnumber_of_isnvs" > {output.isnvs_summary};
        for _file in {input.vcf_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .isnvs.vcf.gz);
            isnv_count=$(bcftools view $_file | grep -v "^#" | wc -l);
            echo -e "$sample\\t$segment\\t$isnv_count" >> {output.isnvs_summary};
        done
        """

rule align_consensus_to_reference_genome:
    input:
        rules.generate_multiqc_report.output.multiqc_report,
        consensus_files = expand(
            rules.rename_sequences.output.consensus_renamed,
            sample=config["samples"],
            allow_missing=True
        )
    output:
        aln_consensus = config['output'] + "assembly/{segment}/consensus/final_consensus/aln.consensus.fasta"
    params:
        path_consensus = config['output'] + "assembly/{segment}/consensus/final_consensus/",
        reference = get_segment_reference
    log:
        config['output'] + "assembly/{segment}/logs/consensus_illumina/align_consensus_to_reference_genome/align_consensus_to_reference_genome.log"
    benchmark:
        config['output'] + "assembly/{segment}/logs/consensus_illumina/align_consensus_to_reference_genome/align_consensus_to_reference_genome.benchmark.txt"
    shell:
        """
        cat {params.reference} {params.path_consensus}/*.renamed.fasta > {params.path_consensus}/consensus.fasta; 
        minimap2 -a --sam-hit-only --secondary=no --score-N=0 {params.reference} {params.path_consensus}/consensus.fasta -o {params.path_consensus}/aln.consensus.sam; 
        gofasta sam toMultiAlign --pad -s {params.path_consensus}/aln.consensus.sam -o {output.aln_consensus}; 
        sed '/^>/ ! s/-/N/g' {output.aln_consensus} > {params.path_consensus}/aln.consensus.indelsMasked.fasta
        """

rule organize_files:
    input:
        fastp_reports = expand(rules.perform_qc.output.html, sample=config["samples"]),
        vcf_files = expand(
            rules.generate_vcf_consensus.output.vcf,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        isn_vcf_files = expand(
            rules.detect_isnv.output.vcf,
            sample=config["samples"], segment=SEGMENTS.keys()
        ) if config.get("run_isnv", False) else [],
        stats_summary = expand(
            rules.calculate_assembly_statistics.output.stats_summary,
            sample=config["samples"], segment=SEGMENTS.keys()
        ),
        consensus_files = expand(
            rules.infer_consensus_sequence.output.consensus,
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
        mkdir -p {params.outdir}samples/
        for sample in {params.samples}; do
            for segment in {params.segments}; do
                mkdir -p {params.outdir}samples/$sample/$segment;
            done
        done
        for _file in {input.fastp_reports}; do
            sample=$(basename $_file _fastp.html | sed 's/^trim.//');
            ln -sf $PWD/$_file {params.outdir}samples/$sample/fastp.html;
        done
        for _file in {input.vcf_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .consensus.vcf.gz);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/$segment/consensus.vcf.gz;
            ln -sf $PWD/$_file.tbi {params.outdir}samples/$sample/$segment/consensus.vcf.gz.tbi;
        done
        for _file in {input.isn_vcf_files} ""; do
            if [ -z "$_file" ]; then continue; fi
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .isnvs.vcf.gz);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/$segment/isnvs.vcf.gz;
            ln -sf $PWD/$_file.tbi {params.outdir}samples/$sample/$segment/isnvs.vcf.gz.tbi;
        done
        for _file in {input.stats_summary}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .stats_summary.csv);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/$segment/stats_summary.csv;
        done
        for _file in {input.consensus_files}; do
            outdir="{params.outdir}"; rel=${{_file#$outdir}}; rel=${{rel#assembly/}};
            segment=$(echo \"$rel\" | cut -d'/' -f1);
            sample=$(basename $_file .consensus.fasta);
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
