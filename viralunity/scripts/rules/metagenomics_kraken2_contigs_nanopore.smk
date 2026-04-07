if run_denovo and run_k2_contigs:
    rule run_kraken2_contigs:
        input:
            fasta = get_final_contigs
        output:
            report = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.report.txt",
            outfile = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/results/{sample}.output.txt"
        threads: config["threads"]
        params:
            database = config["kraken2_database"],
            minimum_hit_group = config.get("minimum_hit_group", 4)
        log:
            config["output"] + "logs/kraken2_contigs/{sample}.log"
        benchmark:
            config["output"] + "logs/kraken2_contigs/{sample}.benchmark.txt"
        conda:
            "../envs/taxonomy.yaml"
        shell:
            """
            set -euo pipefail
            if [ ! -s {input.fasta} ]; then
                echo "WARNING: {input.fasta} empty. Creating dummy Kraken2 contigs outputs." >> {log}
                touch {output.report} {output.outfile}
            else
                kraken2 --db {params.database} --threads {threads} --report-minimizer-data \
                    --minimum-hit-group {params.minimum_hit_group} --report {output.report} \
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
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_taxids.py"

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
        conda:
            "../envs/taxonomy.yaml"
        shell:
            """
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
        conda:
            "../envs/utils.yaml"
        script:
            "../python/summarize_krona_taxa.py"

    rule summarize_taxa_kraken2_contigs_all:
        input:
            expand(
                config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/summary/{sample}.taxa.tsv",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/kraken2_contigs_taxa_summary.tsv"
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
