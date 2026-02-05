rule all:
    input:
        config['output'] + "assembly/consensus/final_consensus/aln.consensus.fasta",
        config['output'] + "isnvs/isnvs_summary.tsv",
        config['output'] + "benchmark.tsv"

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample].split()
    return(reads)

rule perform_qc:
    input:
        get_map_input_fastqs,
    output:
        paired_R1 = config['output'] + "qc/data/trim.p.{sample}_R1.fastq.gz",
        paired_R2 = config['output'] + "qc/data/trim.p.{sample}_R2.fastq.gz",
        unpaired_R1 = config['output'] + "qc/data/trim.u.{sample}_R1.fastq.gz",
        unpaired_R2 = config['output'] + "qc/data/trim.u.{sample}_R2.fastq.gz",
        json = config['output'] + "qc/reports/trim.{sample}_fastp.json",
        html = config['output'] + "qc/reports/trim.{sample}_fastp.html"
    params:
        minimum_length = config["minimum_length"],
        adapters = config["adapters"],
        trim_head = config["trim_head"],
        trim_tail = config["trim_tail"],
        cut_front_mean_quality = config["cut_front_mean_quality"],
        cut_tail_mean_quality = config["cut_tail_mean_quality"],
        cut_right_window_size = config["cut_right_window_size"],
        cut_right_mean_quality = config["cut_right_mean_quality"],
        adapter_args = lambda w: "--detect_adapter_for_pe" if not config["adapters"] else "--adapter_fasta " + config["adapters"]
    threads: config["threads"]
    log:
        config['output'] + "logs/fastp/{sample}.log"
    benchmark:
        config['output'] + "logs/fastp/{sample}.benchmark.txt"
    shell:
        """
        fastp \
            -i {input[0]} \
            -I {input[1]} \
            -o {output.paired_R1} \
            -O {output.paired_R2} \
            --unpaired1 {output.unpaired_R1} \
            --unpaired2 {output.unpaired_R2} \
            --length_required {params.minimum_length} \
            --trim_front1 {params.trim_head} \
            --trim_tail1 {params.trim_tail} \
            --trim_front2 {params.trim_head} \
            --trim_tail2 {params.trim_tail} \
            --cut_front \
            --cut_front_mean_quality {params.cut_front_mean_quality} \
            --cut_tail \
            --cut_tail_mean_quality {params.cut_tail_mean_quality} \
            --cut_right \
            --cut_right_window_size {params.cut_right_window_size} \
            --cut_right_mean_quality {params.cut_right_mean_quality} \
            --json {output.json} \
            --html {output.html} \
            --thread {threads} \
            {params.adapter_args}
        """

rule map_reads:
    input:
        reference = config["reference"],
        R1 = rules.perform_qc.output.paired_R1,
        R2 = rules.perform_qc.output.paired_R2,
    output:
        bam = config['output'] + "assembly/mapped_reads/raw/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/mapped_reads/raw/{sample}.sorted.bam.bai"
    log:
        config['output'] + "logs/minimap2/{sample}.log"
    benchmark:
        config['output'] + "logs/minimap2/{sample}.benchmark.txt"
    threads: config["threads"] 
    shell:
        """
        minimap2 -a -t {threads} -x sr {input.reference} {input.R1} {input.R2} | \
            samtools view -bS -F 4 - | \
            samtools sort -o {output.bam} - 2> {log}
        samtools index {output.bam} {output.bam_index}
        """

rule trim_primer_sequences:
    input:
        bam = rules.map_reads.output.bam,
        bam_index = rules.map_reads.output.bam_index
    output:
        bam = config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam",
        bam_index = config['output'] + "assembly/mapped_reads/trimmed/{sample}.sorted.bam.bai",
        trimmed_info = config['output'] + "assembly/mapped_reads/trimmed/{sample}.trimmed.txt"
    params:
        bed = config["scheme"],
        minimum_length = config["minimum_length"],
        path = config['output'] + "assembly/mapped_reads/raw/"
    log:
        config['output'] + "logs/samtools/ampliconclip/{sample}.log"
    benchmark:
        config['output'] + "logs/samtools/ampliconclip/{sample}.benchmark.txt"
    threads: config["threads"]
    shell:
        """
        if [ {params.bed} == NA ]; then
            mv {input.bam} {output.bam};
            mv {input.bam_index} {output.bam_index};
            touch {output.trimmed_info};
            echo "No primer scheme detected, assuming sequence data came from untargeted sequencing approach (fastq files moved to trimmed dir for convenience)." > {params.path}notes.txt
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

rule detect_isnv:
    input:
        reference = config["reference"],
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        bam = temp(config['output'] + "isnvs/{sample}.lofreq.sorted.bam"),
        bam_index = temp(config['output'] + "isnvs/{sample}.lofreq.sorted.bam.bai"),
        vcf = config['output'] + "isnvs/{sample}.isnvs.vcf.gz",
        vcf_index = config['output'] + "isnvs/{sample}.isnvs.vcf.gz.tbi"
    params:
        af_min_threshold = config["af_isnv_threshold"]
    log:
        config['output'] + "logs/lofreq/{sample}.log"
    benchmark:
        config['output'] + "logs/lofreq/{sample}.benchmark.txt"
    threads: config["threads"]
    shell:
        """
        lofreq indelqual \
            -f {input.reference} \
            -o {output.bam} \
            --dindel \
            {input.bam}
        samtools index {output.bam} {output.bam_index}

        lofreq call-parallel \
            --pp-threads {threads} \
            --call-indels \
            -f {input.reference} \
            -o tmp.vcf \
            {output.bam}

        bcftools view -i 'INFO/AF<0.5 & INFO/AF>={params.af_min_threshold}' tmp.vcf -Oz -o {output.vcf}
        tabix {output.vcf}
        rm tmp.vcf
        """

rule infer_consensus_sequence:
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        consensus = config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.fasta"
    params:
        minimum_depth = config["minimum_depth"],
        af_threshold = config["af_threshold"]
    benchmark:
        config['output'] + "logs/samtools/consensus/{sample}.benchmark.txt"
    log:
        config['output'] + "logs/samtools/consensus/{sample}.log"
    shell:
        "samtools consensus -a -d {params.minimum_depth} -m simple -q -c {params.af_threshold} --show-ins yes {input.bam} -o {output.consensus}"

rule generate_vcf_consensus:
    input:
        reference = config["reference"],
        consensus = rules.infer_consensus_sequence.output.consensus
    output:
        vcf = config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.vcf.gz",
        vcf_index = config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.vcf.gz.tbi"
    benchmark:
        config['output'] + "logs/gsaalign/{sample}.benchmark.txt"
    log:
        config['output'] + "logs/gsaalign/{sample}.log"
    shell:
        """
        out_prefix=$(echo {output.vcf} | sed 's/.vcf.gz//')
        GSAlign \
        -r {input.reference} \
        -q {input.consensus} \
        -o $out_prefix \
        -fmt 1 \
        -sen

        bgzip $out_prefix.vcf
        tabix {output.vcf}
        rm $out_prefix.maf
        """

rule calculate_coverage_basewise:
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        table_cov = config['output'] + "assembly/coverage_stats/{sample}.table_cov_basewise.txt"
    log:
        config['output'] + "logs/bedtools/genomecov/{sample}.log"
    benchmark:
        config['output'] + "logs/bedtools/genomecov/{sample}.benchmark.txt"
    threads: config["threads"]
    shell:
        "bedtools genomecov -d -ibam {input.bam} > {output.table_cov}"

rule rename_sequences:
    input:
        consensus = rules.infer_consensus_sequence.output.consensus
    output:
        consensus_renamed = config['output'] + "assembly/consensus/final_consensus/{sample}.consensus.renamed.fasta"
    log:
        config['output'] + "logs/consensus_illumina/rename_sequences/{sample}.log"
    benchmark:
        config['output'] + "logs/consensus_illumina/rename_sequences/{sample}.benchmark.txt"   
    script:
        "rename_sequences.py"

rule calculate_assembly_statistics:
    input:
        get_map_input_fastqs,
        rules.perform_qc.output.paired_R1,
        rules.trim_primer_sequences.output.bam,
        rules.calculate_coverage_basewise.output.table_cov,
        rules.rename_sequences.output.consensus_renamed
    output:
        stats_summary = config['output'] + "assembly/coverage_stats/{sample}.stats_summary.csv"
    params:
        minimum_depth = config["minimum_depth"]
    log:
        config['output'] + "logs/consensus_illumina/calculate_assembly_statistics/{sample}.log"
    benchmark:
        config['output'] + "logs/consensus_illumina/calculate_assembly_statistics/{sample}.benchmark.txt"   
    script:
        "calculate_assembly_stats.py"

rule unify_assembly_statistics_reports:
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
        echo \"sample_name,number_of_reads,number_of_trim_paired_reads,number_of_mapped_reads,average_depth,percentage_above_10x,percentage_above_100x,percentage_above_1000x,horizontal_coverage\" > {output.unified_stats_summary} ;
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
        vcf_files = expand(rules.detect_isnv.output.vcf, sample=config["samples"])
    output:
        isnvs_summary = config['output'] +  "isnvs/isnvs_summary.tsv"
    log:
        config['output'] + "logs/consensus_illumina/summarize_isnvs/summarize_isnvs.log"
    benchmark:
        config['output'] + "logs/consensus_illumina/summarize_isnvs/summarize_isnvs.benchmark.txt"   
    shell:
        """
        echo -e "sample\tnumber_of_isnvs" > {output.isnvs_summary};
        for _file in {input.vcf_files}; do
            sample=$(basename $_file .isnvs.vcf.gz);
            isnv_count=$(bcftools view $_file | grep -v "^#" | wc -l);
            echo -e "$sample\t$isnv_count" >> {output.isnvs_summary};
        done
        """

rule align_consensus_to_reference_genome:
    input:
        rules.generate_multiqc_report.output.multiqc_report
    output:
        aln_consensus = config['output'] + "assembly/consensus/final_consensus/aln.consensus.fasta"
    params:
        path_consensus = config['output'] + "assembly/consensus/final_consensus/",
        reference = config['reference']
    log:
        config['output'] + "logs/consensus_illumina/align_consensus_to_reference_genome/align_consensus_to_reference_genome.log"
    benchmark:
        config['output'] + "logs/consensus_illumina/align_consensus_to_reference_genome/align_consensus_to_reference_genome.benchmark.txt"   
    shell:
        """
        cat {params.reference} {params.path_consensus}/*.fasta > {params.path_consensus}/consensus.fasta; 
        minimap2 -a --sam-hit-only --secondary=no --score-N=0 {params.reference} {params.path_consensus}/consensus.fasta -o {params.path_consensus}/aln.consensus.sam; 
        gofasta sam toMultiAlign --pad -s {params.path_consensus}/aln.consensus.sam -o {output.aln_consensus}; 
        sed '/^>/! s/-/N/g' {output.aln_consensus} > {params.path_consensus}/aln.consensus.indelsMasked.fasta
        """

rule organize_files:
    input:
        fastp_reports = expand(rules.perform_qc.output.html, sample=config["samples"]),
        vcf_files = expand(rules.generate_vcf_consensus.output.vcf, sample=config["samples"]),
        isn_vcf_files = expand(rules.detect_isnv.output.vcf, sample=config["samples"]),
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
        for _file in {input.isn_vcf_files}; do
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
        
        echo -e "sample\ttask\tseconds\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time" > {output}
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