if run_denovo and run_diamond_contigs:
    rule run_diamond_contigs:
        input:
            fasta = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa",
            db = diamond_db_file
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
        conda:
            "../envs/taxonomy.yaml"
        shell:
            r"""
            set -euo pipefail
            if [ ! -s {input.fasta} ]; then
                echo "WARNING: {input.fasta} empty. Creating dummy DIAMOND contigs output." >> {log}
                touch {output.tsv}
            else
                diamond blastx --db {input.db} --query {input.fasta} \
                    --out {output.tsv} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen qseq qseq_translated full_qseq sseq full_sseq \
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
        conda:
            "../envs/utils.yaml"
        shell:
            """
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
        conda:
            "../envs/utils.yaml"
        shell:
            """
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
        conda:
            "../envs/alignment.yaml"
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
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_diamond_by_idxstats.py"

    rule annotate_diamond_taxonomy:
        input:
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tsv",
            assembly = config["taxids"]
        output:
            annotated = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tax.tsv"
        log:
            config["output"] + "logs/annotate_tax/{sample}.log"
        benchmark:
            config["output"] + "logs/annotate_tax/{sample}.benchmark.txt"
        params:
            taxdump = config["taxdump"]
        conda:
            "../envs/utils.yaml"
        script:
            "../python/annotate_diamond_taxonomy.py"

    rule create_krona_input_from_diamond_contigs:
        input:
            diamond = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.tsv",
            fasta = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa",
            assembly = config["taxids"]
        output:
            krona_input = temp(config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.temp.tsv")
        params:
            data_format = "fasta"
        conda:
            "../envs/utils.yaml"
        script:
            "../python/convert_diamond_output_to_krona_input.py"

    rule filter_krona_input_from_diamond_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.temp.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/results/{sample}.diamond.supported.krona_input.tsv"
        params:
            taxid_column = 1,
            exclude_taxids = EXCLUDE_TAXIDS
        conda:
            "../envs/utils.yaml"
        script:
            "../python/filter_taxids.py"

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
        conda:
            "../envs/taxonomy.yaml"
        shell:
            """
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
        conda:
            "../envs/utils.yaml"
        script:
            "../python/summarize_krona_taxa.py"

    rule summarize_taxa_diamond_contigs_all:
        input:
            expand(
                config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/summary/{sample}.taxa.tsv",
                sample=config["samples"]
            )
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary.tsv"
        conda:
            "../envs/utils.yaml"
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
        conda:
            "../envs/utils.yaml"
        script:
            "../python/add_RPM_to_summary.py"

    rule apply_bleed_filter_diamond_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.bleed.tsv"
        params:
            fraction = config.get("bleed_fraction", 0.005),
            rpm_floor = 1.0,
            rpm_col = "rpm",
        conda:
            "../envs/utils.yaml"
        script:
            "../python/apply_max_rpm_bleed_filter.py"

    rule apply_negative_background_diamond_contigs:
        input:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.bleed.tsv"
        output:
            config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.bleed.neg.tsv"
        params:
            negatives = config.get("negative_controls", []),
            count_col = "mapped_reads",
            p_threshold = config.get("negative_p_threshold", 0.01)
        conda:
            "../envs/utils.yaml"
        script:
            "../python/apply_negative_background_filter.py"
