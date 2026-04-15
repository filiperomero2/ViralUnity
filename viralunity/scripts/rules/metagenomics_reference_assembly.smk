import pandas as pd
import os

def get_checkpoint_inputs(wildcards):
    targets = []
    if config.get("run_kraken2_reads", True):
        targets.append(config["output"] + "metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary.tsv")
    if config.get("run_denovo_assembly", False) and config.get("run_kraken2_contigs", True):
        targets.append(config["output"] + "metagenomics/taxonomic_assignments/kraken2_contigs/kraken2_contigs_taxa_summary.tsv")
    if config.get("run_diamond_reads", False):
        targets.append(config["output"] + "metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary.tsv")
    if config.get("run_denovo_assembly", False) and config.get("run_diamond_contigs", False):
        targets.append(config["output"] + "metagenomics/taxonomic_assignments/diamond_contigs/diamond_contigs_taxa_summary.tsv")
    return targets

checkpoint select_references_meta:
    input:
        get_checkpoint_inputs
    output:
        tsv = config["output"] + "reference_targets.tsv"
    conda:
        "../envs/genome_selection.yaml"
    params:
        summary_dir = config["output"] + "metagenomics/taxonomic_assignments/",
        method = config.get("ref_assembly_method", "kraken2"),
        source = config.get("ref_assembly_source", "reads"),
        reads_count = config.get("ref_assembly_reads_count", 100),
        contigs_count = config.get("ref_assembly_contigs_count", 1),
        families = ",".join(config.get("ref_assembly_families", ["Coronaviridae"])),
        strategy = config.get("ref_selection_strategy", "taxid"),
        genome2taxid = config.get("viral_taxids", "databases/virus_genomes/genome2taxid.tsv"),
        blast_db = config.get("viral_genomes", "databases/virus_genomes/viral.genomes.fasta"),
        blast_qcov = config.get("ref_blast_qcov", 80),
        blast_pident = config.get("ref_blast_pident", 80),
        contigs_dir = config["output"] + "denovo_assembly/viral_contigs/"
    script:
        "../python/select_reference_genomes.py"

# Extract fasta for the given ref_key (family_accession) from the viral genomes database
# ref_key is a unique assembly target identifier: "{family}_{accession}"
rule extract_reference_fasta:
    input:
        tsv = rules.select_references_meta.output.tsv
    output:
        fasta = config["output"] + "assembly/{ref_key}/references/{sample}.fasta"
    params:
        db = config.get("viral_genomes", "databases/virus_genomes/viral.genomes.fasta")
    log:
        config["output"] + "assembly/{ref_key}/logs/extract_reference/{sample}.log"
    conda:
        "../envs/utils.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})
        # Parse the TSV to get the accession for the ref_key/sample
        acc=$(awk -F'\\t' -v s="{wildcards.sample}" -v rk="{wildcards.ref_key}" \
            'NR>1 && $1==s && $2==rk {{print $3}}' {input.tsv} | head -n 1)
        if [ -z "$acc" ]; then
            echo "ERROR: no accession found for sample={wildcards.sample} ref_key={wildcards.ref_key}" | tee -a {log} >&2
            exit 1
        fi
        # Exact-match on the first whitespace-delimited token of the FASTA header
        awk -v seq="$acc" 'BEGIN {{RS=">"; FS="\\n"}} ($1==seq || $1~("^"seq" ")) {{print ">"$0}}' \
            {params.db} > {output.fasta}
        if [ ! -s {output.fasta} ]; then
            echo "ERROR: accession $acc not found in {params.db}" | tee -a {log} >&2
            exit 1
        fi
        echo "Extracted reference $acc for {wildcards.sample}/{wildcards.ref_key}" >> {log}
        """

def get_meta_reference(wildcards):
    return config["output"] + f"assembly/{wildcards.ref_key}/references/{wildcards.sample}.fasta"

REFERENCE = get_meta_reference
SEGMENT_WILDCARD = "{ref_key}/"

if config.get("data") == "illumina":
    include: "alignment_illumina.smk"
    include: "consensus_illumina.smk"
else:
    include: "alignment_nanopore.smk"
    include: "consensus_nanopore.smk"


def get_all_reference_assemblies(wildcards):
    checkpoints.select_references_meta.get()
    tsv_path = config["output"] + "reference_targets.tsv"
    if not os.path.exists(tsv_path):
        return []
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except pd.errors.EmptyDataError:
        return []
    if df.empty:
        return []
    targets = []
    for _, row in df.iterrows():
        sample = row["sample"]
        ref_key = row["ref_key"]
        targets.append(
            config["output"]
            + f"assembly/{ref_key}/consensus/final_consensus/{sample}.consensus.fasta"
        )
    return list(set(targets))


rule collect_reference_assemblies:
    input:
        get_all_reference_assemblies
    output:
        config["output"] + "reference_assembly_done.txt"
    log:
        config["output"] + "logs/collect_reference_assemblies.log"
    benchmark:
        config["output"] + "logs/collect_reference_assemblies.benchmark.txt"
    shell:
        "touch {output} && echo 'Reference assembly complete.' > {log}"

