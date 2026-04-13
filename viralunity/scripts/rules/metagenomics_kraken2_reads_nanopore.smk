if run_k2_reads:
    rule run_kraken2_reads:
        input:
            filtered = rules.remove_host_reads.output.filtered
        output:
            report = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.report.txt",
            outfile = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/results/{sample}.output.txt"
        threads: config.get("run_kraken2_reads_cpus", 2)
        resources:
            mem_mb = config.get("run_kraken2_reads_ram", 4) * 1024
        params:
            database = config["kraken2_database"],
            minimum_hit_group = config.get("minimum_hit_group", 4)
        log:
            config["output"] + "logs/kraken2_reads/{sample}.log"
        benchmark:
            config["output"] + "logs/kraken2_reads/{sample}.benchmark.txt"
        conda:
            "../envs/taxonomy.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.report}) $(dirname {log})
            unc_size=$(gzip -l "{input.filtered}" 2>/dev/null | awk 'NR==2 {{print $1}}' || echo 0)
            if [ -z "$unc_size" ] || [ "$unc_size" = "0" ]; then
                echo "WARNING: {input.filtered} empty. Creating dummy Kraken2 READS outputs." > {log}
                : > {output.report}
                : > {output.outfile}
            else
                kraken2 --db {params.database} --threads {threads} --report-minimizer-data \
                    --minimum-hit-group {params.minimum_hit_group} --report {output.report} \
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
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_taxids.py"

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
        conda:
            "../envs/taxonomy.yaml"
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
        conda:
            "../envs/utils.yaml"
        script:
            "../python/summarize_krona_taxa.py"

    rule summarize_taxa_kraken2_reads_all:
        input:
            expand(
                config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/summary/{sample}.taxa.tsv",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary.tsv"
        conda:
            "../envs/utils.yaml"
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
        conda:
            "../envs/utils.yaml"
        script:
            "../python/add_RPM_to_summary.py"

    rule apply_bleed_filter_kraken2_reads:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.tsv"
        params:
            fraction = config.get("bleed_fraction", 0.005),
            rpm_floor = 1.0,
            rpm_col = "rpm",
        conda:
            "../envs/utils.yaml"
        script:
            "../python/apply_max_rpm_bleed_filter.py"

    rule apply_negative_background_kraken2_reads:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.neg.tsv"
        params:
            negatives = config.get("negative_controls", []),
            count_col = "count",
            p_threshold = config.get("negative_p_threshold", 0.01)
        conda:
            "../envs/utils.yaml"
        script:
            "../python/apply_negative_background_filter.py"
