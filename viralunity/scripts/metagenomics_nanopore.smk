##############################################
# Metagenomics Nanopore – optional dehosting,
# Kraken2/Diamond on reads and (optional) assembly + polishing
##############################################

def get_exclude_taxids():
    exclude = []
    if config.get("remove_human_reads", False):
        exclude.append("9606")
    if config.get("remove_unclassified_reads", False):
        exclude.append("0")
    return exclude

EXCLUDE_TAXIDS = get_exclude_taxids()

def get_sample_to_fastq():
    """Map each sample to its host-filtered FASTQ for read counting (RPM)."""
    return {s: config["output"] + "host_filtered/" + s + ".filtered.fastq.gz" for s in config["samples"]}

def get_map_input_fastqs(wildcards):
    """Return list of read files for this sample (single-end long reads)."""
    paths = config["samples"][wildcards.sample]
    if isinstance(paths, str):
        paths = paths.strip().split()
    return paths

host_filtering_enabled = (config.get("host_reference", "NA") not in ("NA", "", None)) or (config.get("deacon_index", "NA") not in ("NA", "", None))
dehost_with_deacon = config.get("deacon_index", "NA") not in ("NA", "", None)
run_denovo = config.get("run_denovo_assembly", False)
run_polish_racon = config.get("run_polish_racon", False)
run_polish_medaka = config.get("run_polish_medaka", False)
run_k2_reads = config.get("run_kraken2_reads", True)
run_k2_contigs = config.get("run_kraken2_contigs", True)
run_diamond_reads = config.get("run_diamond_reads", False)
run_diamond_contigs = config.get("run_diamond_contigs", False)
has_negative_controls = bool(config.get("negative_controls", []))

def get_final_contigs(wildcards):
    """Path to contigs used for classification (polished, racon, or raw MEGAHIT)."""
    base = config["output"] + "denovo_assembly/megahit/{sample}/"
    s = wildcards.sample
    if run_polish_medaka:
        return base.format(sample=s) + "polished.fasta"
    if run_polish_racon:
        return base.format(sample=s) + "racon.fasta"
    return base.format(sample=s) + "final.contigs.fa"

def get_medaka_assembly_input(wildcards):
    """Assembly input for Medaka (racon output or MEGAHIT)."""
    base = config["output"] + "denovo_assembly/megahit/{sample}/"
    s = wildcards.sample
    if run_polish_racon:
        return base.format(sample=s) + "racon.fasta"
    return base.format(sample=s) + "final.contigs.fa"

##############################################
# Targets
##############################################

def _all_inputs():
    targets = []
    if run_k2_reads:
        if has_negative_controls:
            targets.append(config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.neg.tsv")
        else:
            targets.append(config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.tsv")
    if run_denovo and run_k2_contigs:
        targets.append(config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/kraken2_contigs_taxa_summary.tsv")
    if run_diamond_reads:
        if has_negative_controls:
            targets.append(config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.neg.tsv")
        else:
            targets.append(config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.tsv")
    if run_denovo and run_diamond_contigs:
        if has_negative_controls:
            targets.append(config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.bleed.neg.tsv")
        else:
            targets.append(config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.bleed.tsv")
    if not targets:
        targets = [config["output"] + "metagenomics/metagenomics_summary.txt"]
    return targets

rule all:
    input:
        _all_inputs(),

##############################################
# Dehosting (Nanopore long reads)
# Either Deacon (minimizer index) or minimap2 (host FASTA), or none
##############################################

if dehost_with_deacon:

    rule remove_host_reads:
        input:
            reads = lambda wildcards: config["samples"][wildcards.sample],
            index = config["deacon_index"]
        output:
            filtered = config["output"] + "host_filtered/{sample}.filtered.fastq.gz",
            summary = config["output"] + "logs/remove_host/{sample}.deacon_summary.json",
        threads: config["threads"]
        log:
            config["output"] + "logs/remove_host/{sample}.log"
        benchmark:
            config["output"] + "logs/remove_host/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.filtered})" "$(dirname {log})"
            deacon filter -d -t {threads} -s "{output.summary}" "{input.index}" "{input.reads}" -o "{output.filtered}" 2>> "{log}"
            """

elif host_filtering_enabled:

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
            minimap2 -d {output.index} {input.fasta} 2> {log}
            """

    rule remove_host_reads:
        input:
            reads = lambda wildcards: config["samples"][wildcards.sample],
            index = config["host_reference"] + ".mmi"
        output:
            filtered = config["output"] + "host_filtered/{sample}.filtered.fastq.gz"
        threads: config["threads"]
        log:
            config["output"] + "logs/remove_host/{sample}.log"
        benchmark:
            config["output"] + "logs/remove_host/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.filtered})" "$(dirname {log})"
            minimap2 -t {threads} -ax map-ont "{input.index}" "{input.reads}" 2>> "{log}" \
            | samtools view -@ {threads} -u -f 4 - 2>> "{log}" \
            | samtools fastq -@ {threads} - 2>> "{log}" \
            | gzip -c > "{output.filtered}"
            unc_size="$(gzip -l "{output.filtered}" 2>/dev/null | awk 'NR==2 {{print $2}}' || echo 0)"
            if [[ "${{unc_size}}" == "0" ]] || [[ -z "${{unc_size}}" ]]; then
                echo "WARNING: No reads remained after host filtering for sample {wildcards.sample}." >> "{log}"
            fi
            """

else:

    rule remove_host_reads:
        input:
            reads = lambda wildcards: config["samples"][wildcards.sample],
        output:
            filtered = config["output"] + "host_filtered/{sample}.filtered.fastq.gz"
        log:
            config["output"] + "logs/remove_host/{sample}.log"
        benchmark:
            config["output"] + "logs/remove_host/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.filtered})" "$(dirname {log})"
            if [[ "{input.reads}" == *.gz ]]; then
                cp -f "{input.reads}" "{output.filtered}"
            else
                gzip -c "{input.reads}" > "{output.filtered}"
            fi
            echo "No host genome provided — skipping host read removal." > "{log}"
            """

##############################################
# Kraken2 on reads
##############################################

if run_k2_reads:

    rule run_kraken2_reads:
        input:
            filtered = rules.remove_host_reads.output.filtered
        output:
            report = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.report.txt",
            outfile = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.txt"
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
            unc_size=$(gzip -l "{input.filtered}" 2>/dev/null | awk 'NR==2 {{print $1}}' || echo 0)
            if [ -z "$unc_size" ] || [ "$unc_size" = "0" ]; then
                echo "WARNING: {input.filtered} empty. Creating dummy Kraken2 READS outputs." > {log}
                : > {output.report}
                : > {output.outfile}
            else
                kraken2 --db {params.database} --threads {threads} --report-minimizer-data \
                    --minimum-hit-group 3 --report {output.report} \
                    --output {output.outfile} {input.filtered} 2> {log}
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
                config["output"] + "host_filtered/{sample}.filtered.fastq.gz",
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
            fastq = rules.remove_host_reads.output.filtered,
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
            unc_size=$(gzip -l "{input.fastq}" 2>/dev/null | awk 'NR==2 {{print $1}}' || echo 0)
            if [ -z "$unc_size" ] || [ "$unc_size" = "0" ]; then
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
            fastq = rules.remove_host_reads.output.filtered,
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
                config["output"] + "host_filtered/{sample}.filtered.fastq.gz",
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
# De novo assembly (MEGAHIT, single-end) + optional polishing
##############################################

if run_denovo:

    rule run_megahit:
        input:
            reads = rules.remove_host_reads.output.filtered
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
            if [ ! -s {input.reads} ]; then
                echo "Reads empty; skipping MEGAHIT." > {log}
                mkdir -p $(dirname {output.contigs})
                touch {output.contigs}
            else
                megahit -r {input.reads} -o {params.tempdir} --num-cpu-threads {threads} --tmp-dir "$workdir" >> {log} 2>&1
                mv -T {params.tempdir} {params.outdir}
            fi
            """

    if run_polish_racon:

        rule run_racon:
            input:
                assembly = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa",
                reads = rules.remove_host_reads.output.filtered
            output:
                racon_fasta = config["output"] + "denovo_assembly/megahit/{sample}/racon.fasta"
            threads: config["threads"]
            log:
                config["output"] + "logs/racon/{sample}.log"
            benchmark:
                config["output"] + "logs/racon/{sample}.benchmark.txt"
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {log})
                if [ ! -s {input.assembly} ]; then
                    echo "Assembly empty; skipping Racon." > {log}
                    touch {output.racon_fasta}
                else
                    workdir="${{TMPDIR:-/tmp}}/racon_{wildcards.sample}"
                    mkdir -p "$workdir"
                    minimap2 -t {threads} -x map-ont {input.assembly} {input.reads} > "$workdir/align.paf" 2>> {log}
                    racon -t {threads} {input.reads} "$workdir/align.paf" {input.assembly} > {output.racon_fasta} 2>> {log}
                fi
                """

    if run_polish_medaka:

        rule run_medaka:
            input:
                assembly = get_medaka_assembly_input,
                reads = rules.remove_host_reads.output.filtered
            output:
                polished = config["output"] + "denovo_assembly/megahit/{sample}/polished.fasta",
                bam = config["output"] + "medaka_work/{sample}/calls_to_draft.bam"
            threads: config["threads"]
            log:
                config["output"] + "logs/medaka/{sample}.log"
            benchmark:
                config["output"] + "logs/medaka/{sample}.benchmark.txt"
            params:
                outdir = config["output"] + "medaka_work/{sample}"
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {log})
                if [ ! -s {input.assembly} ]; then
                    echo "No assembly to polish for sample {wildcards.sample}." > {log}
                    mkdir -p $(dirname {output.polished})
                    touch {output.polished}
                    mkdir -p {params.outdir}
                    touch {output.bam}
                else
                    mkdir -p {params.outdir}
                    medaka_consensus \
                        -i {input.reads} \
                        -d {input.assembly} \
                        -o {params.outdir} \
                        -g -r 'N' \
                        -t {threads} &> {log}
                    mv {params.outdir}/consensus.fasta {output.polished}
                fi
                """

##############################################
# Kraken2 on contigs (optional)
##############################################

if run_denovo and run_k2_contigs:

    rule run_kraken2_contigs:
        input:
            fasta = get_final_contigs
        output:
            report = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.report.txt",
            outfile = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.txt"
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
            fasta = get_final_contigs,
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
            contigs = get_final_contigs,
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.tsv"
        output:
            ids = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral.ids.txt",
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

    if run_polish_medaka:

        rule bam_sort_index_idxstats_from_medaka:
            input:
                bam = rules.run_medaka.output.bam,
                viral_ids = rules.extract_viral_contigs.output.ids
            output:
                bam = config["output"] + "mapping/viral/{sample}.viral.bam",
                bai = config["output"] + "mapping/viral/{sample}.viral.bam.bai",
                idxstats = config["output"] + "mapping/viral/{sample}.viral.idxstats.txt",
                idxstats_filtered = config["output"] + "mapping/viral/{sample}.viral.idxstats.filtered.txt"
            threads: config["threads"]
            log:
                config["output"] + "logs/bam_idxstats/{sample}.log"
            benchmark:
                config["output"] + "logs/bam_idxstats/{sample}.benchmark.txt"
            shell:
                r"""
                set -euo pipefail
                mkdir -p "$(dirname {output.bam})" "$(dirname {log})"
                if [ ! -s "{input.bam}" ]; then
                    echo "No Medaka BAM for {wildcards.sample}." > "{log}"
                    touch "{output.bam}" "{output.bai}" "{output.idxstats}" "{output.idxstats_filtered}"
                elif [ ! -s "{input.viral_ids}" ]; then
                    echo "No viral contigs for {wildcards.sample}." > "{log}"
                    touch "{output.bam}" "{output.bai}" "{output.idxstats}" "{output.idxstats_filtered}"
                else
                    refs=()
                    while IFS= read -r r; do [[ -n "$r" ]] && refs+=("$r"); done < "{input.viral_ids}"
                    samtools view -@ {threads} -b "{input.bam}" "$${{refs[@]}}" | samtools sort -@ {threads} -o "{output.bam}" -
                    samtools index -@ {threads} "{output.bam}"
                    samtools idxstats "{output.bam}" > "{output.idxstats}"
                    awk '$3 > 0 && $1 != "*"' "{output.idxstats}" > "{output.idxstats_filtered}"
                fi
                """

    else:

        rule remap_reads_to_viral_contigs:
            input:
                index = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.mmi",
                fasta = config["output"] + "denovo_assembly/viral_contigs/{sample}.viral_contigs.fa",
                reads = rules.remove_host_reads.output.filtered
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
                    minimap2 -t {threads} -x map-ont {input.index} {input.reads} | \
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
            fasta = get_final_contigs,
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
                config["output"] + "host_filtered/{sample}.filtered.fastq.gz",
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
# Fallback summary when no targets enabled
##############################################

rule generate_metagenomics_summary:
    output:
        config["output"] + "metagenomics/metagenomics_summary.txt"
    shell:
        r"""
        mkdir -p $(dirname {output})
        echo "Metagenomics summary (nanopore). Enable run_kraken2_reads, run_diamond_reads, or assembly targets for full outputs." > {output}
        """
