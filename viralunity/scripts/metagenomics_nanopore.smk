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

diamond_db_input_path = config.get("diamond_database", "NA")
if diamond_db_input_path != "NA":
    diamond_db_is_ready = diamond_db_input_path.endswith(".dmnd")
    diamond_db_file = diamond_db_input_path if diamond_db_is_ready else diamond_db_input_path + ".dmnd"
else:
    diamond_db_is_ready = False
    diamond_db_file = "NA"

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

def _all_inputs():
    targets = [
        config["output"] + "benchmark.tsv"
    ]
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

    if config.get("run_reference_assembly", False):
        targets.append(config["output"] + "reference_assembly_done.txt")
    return targets

rule all:
    input:
        _all_inputs(),

if (run_diamond_reads or run_diamond_contigs) and not diamond_db_is_ready and diamond_db_input_path != "NA":
    rule create_diamond_db_shared:
        input:
            diamond_db_input_path
        output:
            diamond_db_file
        log:
            config["output"] + "logs/diamond/diamond_makedb.log"
        benchmark:
            config["output"] + "logs/diamond/diamond_makedb.benchmark.log"
        conda:
            "envs/taxonomy.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output}) $(dirname {log})
            diamond makedb --in {input} --db {output} 2> {log}
            """
include: "rules/metagenomics_dehost_nanopore.smk"
include: "rules/metagenomics_kraken2_reads_nanopore.smk"
include: "rules/metagenomics_diamond_reads_nanopore.smk"
include: "rules/metagenomics_assembly_nanopore.smk"
include: "rules/metagenomics_kraken2_contigs_nanopore.smk"
include: "rules/metagenomics_diamond_contigs_nanopore.smk"
if config.get("run_reference_assembly", False):
    include: "rules/metagenomics_reference_assembly.smk"

rule organize_files:
    conda:
        "envs/utils.yaml"
    input:
        kraken2_reads_reports = expand(rules.run_kraken2_reads.output.report, sample=config["samples"]) if run_k2_reads else [],
        kraken2_reads_krona = expand(rules.create_krona_report_from_kraken2_reads.output, sample=config["samples"]) if run_k2_reads else [],
        diamond_reads_tsv = expand(rules.run_diamond_reads.output.tsv, sample=config["samples"]) if run_diamond_reads else [],
        diamond_reads_krona = expand(rules.create_krona_report_diamond_reads.output, sample=config["samples"]) if run_diamond_reads else [],
        host_filtered = expand(rules.remove_host_reads.output.filtered, sample=config["samples"]),
        megahit_contigs = expand(rules.run_megahit.output.contigs, sample=config["samples"]) if run_denovo else [],
        kraken2_contigs_reports = expand(rules.run_kraken2_contigs.output.report, sample=config["samples"]) if run_denovo and run_k2_contigs else [],
        kraken2_contigs_krona = expand(rules.create_krona_report_from_kraken2_contigs.output, sample=config["samples"]) if run_denovo and run_k2_contigs else [],
        diamond_contigs_tsv = expand(rules.run_diamond_contigs.output.tsv, sample=config["samples"]) if run_denovo and run_diamond_contigs else [],
        diamond_contigs_krona = expand(rules.create_krona_report_diamond_contigs.output, sample=config["samples"]) if run_denovo and run_diamond_contigs else [],
        summary_k2_reads = config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.tsv" if run_k2_reads else [],
        summary_diamond_reads = config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.tsv" if run_diamond_reads else [],
        summary_k2_contigs = config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/kraken2_contigs_taxa_summary.tsv" if run_denovo and run_k2_contigs else [],
        summary_diamond_contigs = config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary_RPM.bleed.tsv" if run_denovo and run_diamond_contigs else [],
    output:
        benchmark = config["output"] + "benchmark.tsv"
    params:
        outdir = config["output"],
        samples = list(config["samples"].keys()),
        run_k2_reads = run_k2_reads,
        run_diamond_reads = run_diamond_reads,
        run_denovo = run_denovo,
        run_k2_contigs = run_k2_contigs,
        run_diamond_contigs = run_diamond_contigs
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}samples/
        for sample in {params.samples}; do
            mkdir -p {params.outdir}samples/$sample;
        done

        # Kraken2 Reads
        if [ "{params.run_k2_reads}" = "True" ]; then
            for _file in {input.kraken2_reads_reports}; do
                sample=$(basename $_file .report.txt);
                ln -sf $PWD/$_file {params.outdir}samples/$sample/kraken2_reads.report.txt;
            done
            for _file in {input.kraken2_reads_krona}; do
                sample=$(basename $_file .output.krona.html);
                ln -sf $PWD/$_file {params.outdir}samples/$sample/kraken2_reads.krona.html;
            done
        fi

        # DIAMOND Reads
        if [ "{params.run_diamond_reads}" = "True" ]; then
            for _file in {input.diamond_reads_tsv}; do
                sample=$(basename $_file .diamond.tsv);
                ln -sf $PWD/$_file {params.outdir}samples/$sample/diamond_reads.tsv;
            done
            for _file in {input.diamond_reads_krona}; do
                sample=$(basename $_file .diamond.krona.html);
                ln -sf $PWD/$_file {params.outdir}samples/$sample/diamond_reads.krona.html;
            done
        fi

        # De novo assembly
        if [ "{params.run_denovo}" = "True" ]; then
            for _file in {input.megahit_contigs}; do
                sample=$(basename $(dirname $_file));
                ln -sf $PWD/$_file {params.outdir}samples/$sample/contigs.fa;
            done
        fi

        # Taxonomy (Contigs)
        if [ "{params.run_denovo}" = "True" ]; then
            if [ "{params.run_k2_contigs}" = "True" ]; then
                for _file in {input.kraken2_contigs_reports}; do
                    sample=$(basename $_file .report.txt);
                    ln -sf $PWD/$_file {params.outdir}samples/$sample/kraken2_contigs.report.txt;
                done
                for _file in {input.kraken2_contigs_krona}; do
                    sample=$(basename $_file .output.krona.html);
                    ln -sf $PWD/$_file {params.outdir}samples/$sample/kraken2_contigs.krona.html;
                done
            fi
            if [ "{params.run_diamond_contigs}" = "True" ]; then
                for _file in {input.diamond_contigs_tsv}; do
                    sample=$(basename $_file .diamond.tsv);
                    ln -sf $PWD/$_file {params.outdir}samples/$sample/diamond_contigs.tsv;
                done
                for _file in {input.diamond_contigs_krona}; do
                    sample=$(basename $_file .diamond.supported.krona.html);
                    ln -sf $PWD/$_file {params.outdir}samples/$sample/diamond_contigs.krona.html;
                done
            fi
        fi

        # Host filtered reads (Nanopore)
        for _file in {input.host_filtered}; do
            sample=$(basename $_file .filtered.fastq.gz);
            ln -sf $PWD/$_file {params.outdir}samples/$sample/host_filtered_reads.fastq.gz;
        done

        # Benchmark aggregation
        echo -e "sample\\ttask\\tseconds\\th:m:s\\tmax_rss\\tmax_vms\\tmax_uss\\tmax_pss\\tio_in\\tio_out\\tmean_load\\tcpu_time" > {output}
        find {params.outdir} -name "*.benchmark.txt" | while read -r file; do
            task=$(basename $(dirname $file))
            sample=$(basename "$file" .benchmark.txt)
            
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
                sample=$(echo "$sample" | sed 's/sample-//')
            fi

            tail -n +2 "$file" | awk -v sample="$sample" -v task="$task" '{{print sample"\\t"task"\\t"$0}}' >> {output}
        done
        """
