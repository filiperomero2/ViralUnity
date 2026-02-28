##############################################
# Metagenomics Illumina – QC, optional dehosting,
# Kraken2/Diamond on reads and (optional) contigs
##############################################

def get_exclude_taxids():
    """TaxIDs to exclude from classification outputs (e.g. human 9606, unclassified 0)."""
    exclude = []
    if config.get("remove_human_reads", False):
        exclude.append("9606")
    if config.get("remove_unclassified_reads", False):
        exclude.append("0")
    return exclude

EXCLUDE_TAXIDS = get_exclude_taxids()

def get_sample_to_fastq():
    """Map each sample to its merged host-filtered FASTQ for read counting (RPM)."""
    return {s: config["output"] + "host_filtered/" + s + ".merged.fastq.gz" for s in config["samples"]}

def get_map_input_fastqs(wildcards):
    """Return [R1, R2] for Trimmomatic / host filtering (Illumina paired-end)."""
    paths = config["samples"][wildcards.sample].strip().split()
    return paths

def get_final_input_fastq(wildcards):
    """Single merged FASTQ per sample for read-level Kraken2/Diamond."""
    return config["output"] + "host_filtered/" + wildcards.sample + ".merged.fastq.gz"

host_filtering_enabled = config.get("host_reference", "NA") != "NA"
run_denovo = config.get("run_denovo_assembly", False)
run_k2_reads = config.get("run_kraken2_reads", True)
run_k2_contigs = config.get("run_kraken2_contigs", True)
run_diamond_reads = config.get("run_diamond_reads", False)
run_diamond_contigs = config.get("run_diamond_contigs", False)

##############################################
# Targets
##############################################

def _all_inputs():
    targets = [config["output"] + "qc/reports/multiqc_report.html"]
    if run_k2_reads:
        targets.append(config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.tsv")
    if run_denovo and run_k2_contigs:
        targets.append(config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/kraken2_contigs_taxa_summary.tsv")
    if run_diamond_reads:
        targets.append(config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.tsv")
    if run_denovo and run_diamond_contigs:
        targets.append(config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.bleed.tsv")
    return targets

rule all:
    input:
        _all_inputs(),

##############################################
# QC and trimming (Illumina) – fastp
##############################################

def get_fastp_adapter_args(wildcards):
    """Use --detect_adapter_for_pe when no adapters file is provided."""
    adapters = config.get("adapters", "NA")
    if not adapters or str(adapters).strip() in ("", "NA"):
        return "--detect_adapter_for_pe"
    return "--adapter_fasta " + str(adapters).strip()

rule perform_qc:
    input:
        get_map_input_fastqs,
    output:
        paired_R1 = config["output"] + "qc/data/trim.p.{sample}_R1.fastq.gz",
        paired_R2 = config["output"] + "qc/data/trim.p.{sample}_R2.fastq.gz",
        unpaired_R1 = config["output"] + "qc/data/trim.u.{sample}_R1.fastq.gz",
        unpaired_R2 = config["output"] + "qc/data/trim.u.{sample}_R2.fastq.gz",
        json = config["output"] + "qc/reports/trim.{sample}_fastp.json",
        html = config["output"] + "qc/reports/trim.{sample}_fastp.html"
    params:
        minimum_length = config.get("minimum_length", 50),
        trim_head = config.get("trim_head", 0),
        trim_tail = config.get("trim_tail", 0),
        cut_front_mean_quality = config.get("cut_front_mean_quality", 20),
        cut_tail_mean_quality = config.get("cut_tail_mean_quality", 20),
        cut_right_window_size = config.get("cut_right_window_size", 4),
        cut_right_mean_quality = config.get("cut_right_mean_quality", 20),
        adapter_args = get_fastp_adapter_args,
    threads: config["threads"]
    log:
        config["output"] + "logs/fastp/{sample}.log"
    benchmark:
        config["output"] + "logs/fastp/{sample}.benchmark.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.paired_R1}) $(dirname {output.json}) $(dirname {log})
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
            {params.adapter_args} 2> {log}
        """

rule generate_qc_report:
    input:
        config["output"] + "qc/data/trim.p.{sample}_R1.fastq.gz",
        config["output"] + "qc/data/trim.p.{sample}_R2.fastq.gz"
    output:
        config["output"] + "qc/reports/trim.p.{sample}_R1_fastqc.html",
        config["output"] + "qc/reports/trim.p.{sample}_R2_fastqc.html",
        config["output"] + "qc/reports/trim.p.{sample}_R1_fastqc.zip",
        config["output"] + "qc/reports/trim.p.{sample}_R2_fastqc.zip"
    threads: config["threads"]
    params:
        temp = config["output"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.temp}/qc/reports
        fastqc -q -t {threads} -o {params.temp}/qc/reports/ {input}
        """

##############################################
# Dehosting (Illumina paired-end)
##############################################

if host_filtering_enabled:

    rule index_host_genome:
        input:
            fasta = config["host_reference"]
        output:
            index = config["host_reference"] + ".mmi"
        threads: config["threads"]
        log:
            config["output"] + "logs/remove_host/host_genome_indexing.log"
        benchmark:
            config["output"] + "logs/remove_host/host_genome_indexing.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            minimap2 -d {output.index} {input.fasta} > {log} 2>&1
            """

    rule remove_host_reads:
        input:
            paired_R1 = rules.perform_qc.output.paired_R1,
            paired_R2 = rules.perform_qc.output.paired_R2,
            index = config["host_reference"] + ".mmi"
        output:
            filtered_R1 = config["output"] + "host_filtered/{sample}.R1.filtered.fastq.gz",
            filtered_R2 = config["output"] + "host_filtered/{sample}.R2.filtered.fastq.gz",
        threads: config["threads"]
        log:
            config["output"] + "logs/remove_host/{sample}.log"
        benchmark:
            config["output"] + "logs/remove_host/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.filtered_R1}) $(dirname {log})
            tmp_bam=$(mktemp --suffix .unmapped.bam)
            minimap2 -t {threads} -ax sr {input.index} {input.paired_R1} {input.paired_R2} \
            | samtools view -@ {threads} -b -f 4 - > "$tmp_bam"
            unmapped=$(samtools view -c "$tmp_bam")
            if [ "$unmapped" -eq 0 ]; then
                touch {output.filtered_R1} {output.filtered_R2}
                echo "All reads mapped to host; empty FASTQs written." >> {log}
            else
                samtools fastq -@ {threads} -1 {output.filtered_R1} -2 {output.filtered_R2} "$tmp_bam" >> {log} 2>&1
            fi
            rm -f "$tmp_bam"
            """

else:

    rule remove_host_reads:
        input:
            paired_R1 = rules.perform_qc.output.paired_R1,
            paired_R2 = rules.perform_qc.output.paired_R2,
        output:
            filtered_R1 = config["output"] + "host_filtered/{sample}.R1.filtered.fastq.gz",
            filtered_R2 = config["output"] + "host_filtered/{sample}.R2.filtered.fastq.gz",
        log:
            config["output"] + "logs/remove_host/{sample}.log"
        benchmark:
            config["output"] + "logs/remove_host/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.filtered_R1}) $(dirname {log})
            cp {input.paired_R1} {output.filtered_R1}
            cp {input.paired_R2} {output.filtered_R2}
            echo "No host genome provided — skipping host read removal." > {log}
            """

##############################################
# Merge host-filtered R1+R2 for read-level tools
##############################################

rule merge_host_filtered_reads:
    input:
        filtered_R1 = rules.remove_host_reads.output.filtered_R1,
        filtered_R2 = rules.remove_host_reads.output.filtered_R2,
    output:
        merged = temp(config["output"] + "host_filtered/{sample}.merged.fastq.gz"),
    log:
        config["output"] + "logs/merge_reads/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.merged}) $(dirname {log})
        decompress() {{
            case "$1" in
                *.gz) gzip -dc "$1" 2>/dev/null || true ;;
                *) [ -s "$1" ] && cat "$1" || true ;;
            esac
        }}
        {{ decompress "{input.filtered_R1}"; decompress "{input.filtered_R2}"; }} | gzip -c > {output.merged}
        """

##############################################
# Kraken2 on reads
##############################################

rule run_kraken2_reads:
    input:
        filtered_R1 = rules.remove_host_reads.output.filtered_R1,
        filtered_R2 = rules.remove_host_reads.output.filtered_R2,
    output:
        report = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.report.txt",
        outfile = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.txt",
    threads: config["threads"]
    params:
        database = config["kraken2_database"]
    log:
        config["output"] + "logs/kraken2_reads/{sample}.log"
    benchmark:
        config["output"] + "logs/kraken2_reads/{sample}.benchmark.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.report}) $(dirname {log})
        unc_size_R1=$(gzip -l "{input.filtered_R1}" | awk 'NR==2 {{print $1}}')
        unc_size_R2=$(gzip -l "{input.filtered_R2}" | awk 'NR==2 {{print $1}}')
        if [ "$unc_size_R1" = "0" ] || [ "$unc_size_R2" = "0" ]; then
            echo "WARNING: {input.filtered_R1} or {input.filtered_R2} empty. Creating dummy Kraken2 READS outputs." > {log}
            : > {output.report}
            : > {output.outfile}
        else
            kraken2 --db {params.database} --threads {threads} --report-minimizer-data \
                --minimum-hit-group 3 --report {output.report} \
                --output {output.outfile} {input.filtered_R1} {input.filtered_R2} 2> {log}
        fi
        """

rule create_krona_input_from_kraken2_reads:
    input:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.txt"
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.krona.txt"
    params:
        keep_columns = [1, 2],
        taxid_column = 1,
        exclude_taxids = EXCLUDE_TAXIDS
    script:
        "filter_taxids.py"

rule create_krona_report_from_kraken2_reads:
    input:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.krona.txt"
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/reports/{sample}.output.krona.html"
    params:
        krona_database = config["krona_database"]
    log:
        config["output"] + "logs/krona_kraken2_reads/{sample}.log"
    benchmark:
        config["output"] + "logs/krona_kraken2_reads/{sample}.benchmark.txt"
    shell:
        r"""
        set -euo pipefail
        if [ -s {input} ]; then
            ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}
        else
            echo "Empty krona input (kraken2 reads)." >> {log}
            touch {output}
        fi
        """

rule summarize_taxa_kraken2_reads:
    input:
        krona = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.krona.txt",
        plot = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/reports/{sample}.output.krona.html"
    output:
        temp(config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/summary/{sample}.taxa.tsv")
    params:
        taxdump = config["taxdump"],
        tool = "kraken2",
        mode = "reads",
        sample = "{sample}"
    script:
        "summarize_krona_taxa.py"

rule summarize_taxa_kraken2_reads_all:
    input:
        expand(
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/summary/{sample}.taxa.tsv",
            sample=config["samples"]
        )
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary.tsv"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output})
        header_file=""
        for f in {input}; do [ -s "$f" ] && header_file="$f" && break; done
        if [ -n "$header_file" ]; then head -n 1 "$header_file" > {output}; else echo -e "sample\ttool\tmode\trank\ttaxid\tname\tcount\tpercent\tsource" > {output}; fi
        for f in {input}; do [ -s "$f" ] && tail -n +2 "$f" >> {output}; done
        """

rule add_RPM_to_kraken2_reads_summary:
    input:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary.tsv",
        merged_fastqs = expand(
            config["output"] + "host_filtered/{sample}.merged.fastq.gz",
            sample=config["samples"]
        )
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.tsv"
    params:
        sample_to_fastq = get_sample_to_fastq(),
        reads_col = "count"
    script:
        "add_RPM_to_summary.py"

rule apply_bleed_filter_kraken2_reads:
    input:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.tsv"
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.tsv"
    params:
        fraction = config.get("bleed_fraction", 0.005),
        rpm_floor = 1.0,
        rpm_col = "rpm",
    script:
        "apply_max_rpm_bleed_filter.py"

rule apply_negative_background_kraken2_reads:
    input:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.tsv"
    output:
        config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.neg.tsv"
    params:
        negatives = config.get("negative_controls", []),
        count_col = "count",
        p_threshold = config.get("negative_p_threshold", 0.01)
    script:
        "apply_negative_background_filter.py"

##############################################
# DIAMOND on reads (optional)
##############################################

if run_diamond_reads:

    rule create_diamond_db:
        input:
            config["diamond_database"]
        output:
            config["diamond_database"] + ".dmnd"
        log:
            config["output"] + "logs/diamond/diamond_makedb.log"
        benchmark:
            config["output"] + "logs/diamond/diamond_makedb.benchmark.log"
        shell:
            r"""
            set -euo pipefail
            if [ "{input}" = "NA" ] || [ -z "{input}" ]; then
                echo "DIAMOND database not provided." > {log}
                exit 1
            fi
            diamond makedb --in {input} --db {output} 2> {log}
            """

    rule run_diamond_reads:
        input:
            fastq = rules.merge_host_filtered_reads.output.merged,
            db = config["diamond_database"] + ".dmnd"
        output:
            tsv = config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.tsv"
        threads: config["threads"]
        log:
            config["output"] + "logs/diamond_reads/{sample}.log"
        benchmark:
            config["output"] + "logs/diamond_reads/{sample}.benchmark.txt"
        params:
            sensitivity = config.get("diamond_sensitivity", "sensitive"),
            evalue = config.get("evalue", 0.001)
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})
            unc_size=$(gzip -l "{input.fastq}" | awk 'NR==2 {{print $1}}')
            if [ "$unc_size" = "0" ]; then
                echo "WARNING: {input.fastq} empty. Creating dummy DIAMOND READS output." > {log}
                touch {output.tsv}
            else
                diamond blastx --db {input.db} --query {input.fastq} \
                    --out {output.tsv} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
                    --max-target-seqs 1 --evalue {params.evalue} \
                    --{params.sensitivity} --threads {threads} 2> {log}
            fi
            """

    rule create_krona_input_from_diamond_reads:
        input:
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.tsv",
            fastq = rules.merge_host_filtered_reads.output.merged,
            assembly = config["assembly_summary"]
        output:
            krona_input = temp(config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.krona_input.temp.tsv")
        params:
            data_format = "fastq"
        script:
            "convert_diamond_output_to_krona_input.py"

    rule filter_krona_input_from_diamond_reads:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.krona_input.temp.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.supported.krona_input.tsv"
        params:
            taxid_column = 1,
            exclude_taxids = EXCLUDE_TAXIDS
        script:
            "filter_taxids.py"

    rule create_krona_report_diamond_reads:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.supported.krona_input.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/reports/{sample}.diamond.krona.html"
        params:
            krona_database = config["krona_database"]
        log:
            config["output"] + "logs/krona_diamond_reads/{sample}.log"
        benchmark:
            config["output"] + "logs/krona_diamond_reads/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            if [ -s {input} ]; then
                ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}
            else
                echo "Empty krona input (diamond reads)." >> {log}
                touch {output}
            fi
            """

    rule summarize_taxa_diamond_reads:
        input:
            krona = config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.supported.krona_input.tsv",
            plot = config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/reports/{sample}.diamond.krona.html"
        output:
            temp(config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/summary/{sample}.taxa.tsv")
        params:
            taxdump = config["taxdump"],
            tool = "diamond",
            mode = "reads",
            sample = "{sample}"
        script:
            "summarize_krona_taxa.py"

    rule summarize_taxa_diamond_reads_all:
        input:
            expand(
                config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/summary/{sample}.taxa.tsv",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary.tsv"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output})
            header_file=""
            for f in {input}; do [ -s "$f" ] && header_file="$f" && break; done
            if [ -n "$header_file" ]; then head -n 1 "$header_file" > {output}; else echo -e "sample\ttool\tmode\trank\ttaxid\tname\tcount\tpercent\tsource" > {output}; fi
            for f in {input}; do [ -s "$f" ] && tail -n +2 "$f" >> {output}; done
            """

    rule add_RPM_to_diamond_reads_summary:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary.tsv",
            merged_fastqs = expand(
                config["output"] + "host_filtered/{sample}.merged.fastq.gz",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.tsv"
        params:
            sample_to_fastq = get_sample_to_fastq(),
            reads_col = "count"
        script:
            "add_RPM_to_summary.py"

    rule apply_bleed_filter_diamond_reads:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.tsv"
        params:
            fraction = config.get("bleed_fraction", 0.005),
            rpm_floor = 1.0,
            rpm_col = "rpm",
        script:
            "apply_max_rpm_bleed_filter.py"

    rule apply_negative_background_diamond_reads:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.neg.tsv"
        params:
            negatives = config.get("negative_controls", []),
            count_col = "count",
            p_threshold = config.get("negative_p_threshold", 0.01)
        script:
            "apply_negative_background_filter.py"

##############################################
# De novo assembly (MEGAHIT)
##############################################

if run_denovo:

    rule run_megahit:
        input:
            filtered_R1 = rules.remove_host_reads.output.filtered_R1,
            filtered_R2 = rules.remove_host_reads.output.filtered_R2,
        output:
            contigs = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa"
        threads: config["threads"]
        log:
            config["output"] + "logs/megahit/{sample}.log"
        benchmark:
            config["output"] + "logs/megahit/{sample}.benchmark.txt"
        params:
            tempdir = temp(config["output"] + "denovo_assembly/megahit/temp_{sample}"),
            outdir = config["output"] + "denovo_assembly/megahit/{sample}"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {log})
            workdir="${{TMPDIR:-/tmp}}/megahit_{wildcards.sample}"
            mkdir -p "$workdir"
            size_R1=$([ -s {input.filtered_R1} ] && echo 1 || echo 0)
            size_R2=$([ -s {input.filtered_R2} ] && echo 1 || echo 0)
            if [ "$size_R1" -eq 1 ] && [ "$size_R2" -eq 1 ]; then
                megahit -1 {input.filtered_R1} -2 {input.filtered_R2} -o {params.tempdir} \
                    --num-cpu-threads {threads} --tmp-dir "$workdir" >> {log} 2>&1
                mv -T {params.tempdir} {params.outdir}
            else
                echo "R1 and/or R2 empty; skipping MEGAHIT." > {log}
                mkdir -p $(dirname {output.contigs})
                touch {output.contigs}
            fi
            """

##############################################
# Kraken2 on contigs (optional)
##############################################

if run_denovo and run_k2_contigs:

    rule run_kraken2_contigs:
        input:
            fasta = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa"
        output:
            report = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.report.txt",
            outfile = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.txt",
        threads: config["threads"]
        params:
            database = config["kraken2_database"]
        log:
            config["output"] + "logs/kraken2_contigs/{sample}.log"
        benchmark:
            config["output"] + "logs/kraken2_contigs/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            if [ ! -s {input.fasta} ]; then
                echo "WARNING: {input.fasta} empty. Creating dummy Kraken2 contigs outputs." >> {log}
                touch {output.report} {output.outfile}
            else
                kraken2 --db {params.database} --threads {threads} --report-minimizer-data \
                    --minimum-hit-group 3 --report {output.report} \
                    --output {output.outfile} {input.fasta} 2> {log}
            fi
            """

    rule create_krona_input_from_kraken2_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.txt"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.krona.txt"
        params:
            keep_columns = [1, 2],
            taxid_column = 1,
            exclude_taxids = EXCLUDE_TAXIDS
        script:
            "filter_taxids.py"

    rule create_krona_report_from_kraken2_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.krona.txt"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/reports/{sample}.output.krona.html"
        params:
            krona_database = config["krona_database"]
        log:
            config["output"] + "logs/krona_kraken2_contigs/{sample}.log"
        benchmark:
            config["output"] + "logs/krona_kraken2_contigs/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            if [ -s {input} ]; then
                ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}
            else
                echo "Empty krona input (kraken2 contigs)." >> {log}
                touch {output}
            fi
            """

    rule summarize_taxa_kraken2_contigs:
        input:
            krona = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.krona.txt",
            plot = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/reports/{sample}.output.krona.html"
        output:
            temp(config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summary/{sample}.taxa.tsv")
        params:
            taxdump = config["taxdump"],
            tool = "kraken2",
            mode = "contigs",
            sample = "{sample}"
        script:
            "summarize_krona_taxa.py"

    rule summarize_taxa_kraken2_contigs_all:
        input:
            expand(
                config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summary/{sample}.taxa.tsv",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/kraken2_contigs_taxa_summary.tsv"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output})
            header_file=""
            for f in {input}; do [ -s "$f" ] && header_file="$f" && break; done
            if [ -n "$header_file" ]; then head -n 1 "$header_file" > {output}; else echo -e "sample\ttool\tmode\trank\ttaxid\tname\tcount\tpercent\tsource" > {output}; fi
            for f in {input}; do [ -s "$f" ] && tail -n +2 "$f" >> {output}; done
            """

##############################################
# DIAMOND on contigs + read support (optional)
##############################################

if run_denovo and run_diamond_contigs:

    rule run_diamond_contigs:
        input:
            fasta = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa",
            db = config["diamond_database"] + ".dmnd"
        output:
            tsv = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.tsv"
        threads: config["threads"]
        log:
            config["output"] + "logs/diamond_contigs/{sample}.log"
        benchmark:
            config["output"] + "logs/diamond_contigs/{sample}.benchmark.txt"
        params:
            sensitivity = config.get("diamond_sensitivity", "sensitive"),
            evalue = config.get("evalue", 0.001)
        shell:
            r"""
            set -euo pipefail
            if [ ! -s {input.fasta} ]; then
                echo "WARNING: {input.fasta} empty. Creating dummy DIAMOND contigs output." >> {log}
                touch {output.tsv}
            else
                diamond blastx --db {input.db} --query {input.fasta} \
                    --out {output.tsv} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq qseq_translated full_qseq sseq full_sseq \
                    --max-target-seqs 1 --evalue {params.evalue} \
                    --{params.sensitivity} --threads {threads} 2> {log}
            fi
            """

    rule extract_viral_contigs:
        input:
            contigs = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa",
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.tsv"
        output:
            ids = temp(config["output"] + "denovo_assembly/viral_contigs/{sample}.viral.ids.txt"),
            fasta = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.fa"
        threads: 1
        log:
            config["output"] + "logs/viral_contigs/{sample}.log"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.fasta}) $(dirname {log})
            cut -f1 {input.diamond} | sort -u > {output.ids}
            if [ -s {output.ids} ]; then
                seqtk subseq {input.contigs} {output.ids} > {output.fasta}
            else
                echo "No viral contigs for {wildcards.sample}" >> {log}
                touch {output.fasta}
            fi
            """

    rule index_viral_contigs:
        input:
            fasta = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.fa"
        output:
            index = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.mmi"
        threads: config["threads"]
        shell:
            r"""
            if [ -s {input.fasta} ]; then
                minimap2 -d {output.index} {input.fasta}
            else
                touch {output.index}
            fi
            """

    rule remap_reads_to_viral_contigs:
        input:
            index = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.mmi",
            fasta = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.fa",
            R1 = rules.remove_host_reads.output.filtered_R1,
            R2 = rules.remove_host_reads.output.filtered_R2
        output:
            bam = config["output"] + "mapping/viral/{sample}.viral.bam",
            bai = config["output"] + "mapping/viral/{sample}.viral.bam.bai",
            idxstats = config["output"] + "mapping/viral/{sample}.viral.idxstats.txt"
        threads: config["threads"]
        log:
            config["output"] + "logs/remap_viral/{sample}.log"
        benchmark:
            config["output"] + "logs/remap_viral/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.bam}) $(dirname {log})
            if [ ! -s {input.fasta} ]; then
                echo "No viral contigs for {wildcards.sample}; skipping remap." >> {log}
                touch {output.bam} {output.bai} {output.idxstats}
            else
                minimap2 -t {threads} -ax sr {input.index} {input.R1} {input.R2} | \
                    samtools sort -@ {threads} -o {output.bam} -
                samtools index -@ {threads} {output.bam}
                samtools idxstats {output.bam} > {output.idxstats}
            fi
            """

    rule diamond_filter_by_idxstats:
        input:
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.tsv",
            idxstats = config["output"] + "mapping/viral/{sample}.viral.idxstats.txt"
        output:
            filtered = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tsv"
        params:
            min_mapped = 1
        log:
            config["output"] + "logs/diamond_filter/{sample}.log"
        benchmark:
            config["output"] + "logs/diamond_filter/{sample}.benchmark.txt"
        script:
            "filter_diamond_by_idxstats.py"

    rule annotate_diamond_taxonomy:
        input:
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tsv",
            assembly = config["assembly_summary"]
        output:
            annotated = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tax.tsv"
        log:
            config["output"] + "logs/annotate_tax/{sample}.log"
        benchmark:
            config["output"] + "logs/annotate_tax/{sample}.benchmark.txt"
        params:
            taxdump = config["taxdump"]
        script:
            "annotate_diamond_taxonomy.py"

    rule create_krona_input_from_diamond_contigs:
        input:
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tsv",
            fasta = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa",
            assembly = config["assembly_summary"]
        output:
            krona_input = temp(config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.temp.tsv")
        params:
            data_format = "fasta"
        script:
            "convert_diamond_output_to_krona_input.py"

    rule filter_krona_input_from_diamond_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.temp.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.tsv"
        params:
            taxid_column = 1,
            exclude_taxids = EXCLUDE_TAXIDS
        script:
            "filter_taxids.py"

    rule create_krona_report_diamond_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/reports/{sample}.diamond.supported.krona.html"
        params:
            krona_database = config["krona_database"]
        log:
            config["output"] + "logs/krona_diamond_contigs/{sample}.log"
        benchmark:
            config["output"] + "logs/krona_diamond_contigs/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            if [ -s {input} ]; then
                ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}
            else
                echo "Empty krona input (diamond contigs)." >> {log}
                touch {output}
            fi
            """

    rule summarize_taxa_diamond_contigs:
        input:
            krona = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.tsv",
            plot = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/reports/{sample}.diamond.supported.krona.html",
            annotated = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tax.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/summary/{sample}.taxa.tsv"
        params:
            taxdump = config["taxdump"],
            tool = "diamond",
            mode = "contigs",
            sample = "{sample}"
        script:
            "summarize_krona_taxa.py"

    rule summarize_taxa_diamond_contigs_all:
        input:
            expand(
                config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/summary/{sample}.taxa.tsv",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary.tsv"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output})
            header_file=""
            for f in {input}; do [ -s "$f" ] && header_file="$f" && break; done
            if [ -n "$header_file" ]; then head -n 1 "$header_file" > {output}; else echo -e "sample\ttool\tmode\trank\ttaxid\tname\tcount\tpercent\tsource\tmapped_reads" > {output}; fi
            for f in {input}; do [ -s "$f" ] && tail -n +2 "$f" >> {output}; done
            """

    rule add_RPM_to_diamond_contigs_summary:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary.tsv",
            merged_fastqs = expand(
                config["output"] + "host_filtered/{sample}.merged.fastq.gz",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.tsv"
        params:
            sample_to_fastq = get_sample_to_fastq(),
            reads_col = "mapped_reads"
        script:
            "add_RPM_to_summary.py"

    rule apply_bleed_filter_diamond_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.bleed.tsv"
        params:
            fraction = config.get("bleed_fraction", 0.005),
            rpm_floor = 1.0,
            rpm_col = "rpm",
        script:
            "apply_max_rpm_bleed_filter.py"

    rule apply_negative_background_diamond_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.bleed.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.bleed.neg.tsv"
        params:
            negatives = config.get("negative_controls", []),
            count_col = "mapped_reads",
            p_threshold = config.get("negative_p_threshold", 0.01)
        script:
            "apply_negative_background_filter.py"

##############################################
# MultiQC (depends on QC reports)
##############################################

rule generate_multiqc_report:
    input:
        expand(
            config["output"] + "qc/reports/trim.p.{sample}_R1_fastqc.html",
            sample=config["samples"]
        )
    output:
        config["output"] + "qc/reports/multiqc_report.html"
    params:
        temp = config["output"]
    shell:
        r"""
        multiqc -f -s -o {params.temp}/qc/reports/ {params.temp}/qc/reports/
        """
