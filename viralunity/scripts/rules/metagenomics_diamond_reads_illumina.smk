if run_diamond_reads:
    rule run_diamond_reads:
        input:
            fastq = rules.merge_host_filtered_reads.output.merged,
            db = diamond_db_file
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
        conda:
            "../envs/taxonomy.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})
            unc_size=$(gzip -l "{input.fastq}" | awk 'NR==2 {{print $1}}')
            if [ "$unc_size" = "0" ]; then
                echo "WARNING: {input.fastq} empty. Creating dummy DIAMOND READS output." > {log}
                touch {output.tsv}
            else
                diamond blastx --db {input.db} --query {input.fastq} \
                    --out {output.tsv} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
                    --max-target-seqs 1 --evalue {params.evalue} \
                    --{params.sensitivity} --threads {threads} 2> {log}
            fi
            """

    rule create_krona_input_from_diamond_reads:
        input:
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.tsv",
            fastq = rules.merge_host_filtered_reads.output.merged,
            assembly = config["taxids"]
        output:
            krona_input = temp(config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.krona_input.temp.tsv")
        params:
            data_format = "fastq"
        conda:
            "../envs/utils.yaml"
        script:
            "../python/convert_diamond_output_to_krona_input.py"

    rule filter_krona_input_from_diamond_reads:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.krona_input.temp.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/results/{sample}.diamond.supported.krona_input.tsv"
        params:
            taxid_column = 1,
            exclude_taxids = EXCLUDE_TAXIDS
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_taxids.py"

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
        conda:
            "../envs/taxonomy.yaml"
        shell:
            """
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
        conda:
            "../envs/utils.yaml"
        script:
            "../python/summarize_krona_taxa.py"

    rule summarize_taxa_diamond_reads_all:
        input:
            expand(
                config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/summary/{sample}.taxa.tsv",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary.tsv"
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
        conda:
            "../envs/utils.yaml"
        script:
            "../python/add_RPM_to_summary.py"

    rule apply_bleed_filter_diamond_reads:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.tsv"
        params:
            fraction = config.get("bleed_fraction", 0.005),
            rpm_floor = 1.0,
            rpm_col = "rpm",
        conda:
            "../envs/utils.yaml"
        script:
            "../python/apply_max_rpm_bleed_filter.py"

    rule apply_negative_background_diamond_reads:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.neg.tsv"
        params:
            negatives = config.get("negative_controls", []),
            count_col = "count",
            p_threshold = config.get("negative_p_threshold", 0.01)
        conda:
            "../envs/utils.yaml"
        script:
            "../python/apply_negative_background_filter.py"
