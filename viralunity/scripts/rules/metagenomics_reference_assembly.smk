import pandas as pd

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
        method = config.get("reference assembly", {}).get("ref_assembly_method", "kraken2"),
        source = config.get("reference assembly", {}).get("ref_assembly_source", "reads"),
        reads_count = config.get("reference assembly", {}).get("ref_assembly_reads_count", 100),
        contigs_count = config.get("reference assembly", {}).get("ref_assembly_contigs_count", 1),
        families = ",".join(config.get("reference assembly", {}).get("ref_assembly_families", ["Coronaviridae"])),
        strategy = config.get("reference assembly", {}).get("ref_selection_strategy", "taxid"),
        genome2taxid = config.get("databases", {}).get("viral_taxids", "databases/virus_genomes/genome2taxid.tsv"),
        blast_db = config.get("databases", {}).get("viral_genomes", "databases/virus_genomes/viral.genomes.fasta"),
        blast_qcov = config.get("reference assembly", {}).get("ref_blast_qcov", 80),
        blast_pident = config.get("reference assembly", {}).get("ref_blast_pident", 80),
        contigs_dir = config["output"] + "denovo_assembly/megahit/"
    shell:
        """
        python scripts/python/select_reference_genomes.py \
            --summary-dir {params.summary_dir} \
            --method {params.method} \
            --source {params.source} \
            --reads-count {params.reads_count} \
            --contigs-count {params.contigs_count} \
            --families '{params.families}' \
            --strategy {params.strategy} \
            --genome2taxid {params.genome2taxid} \
            --blast-db {params.blast_db} \
            --blast-qcov {params.blast_qcov} \
            --blast-pident {params.blast_pident} \
            --contigs-dir {params.contigs_dir} \
            --out-tsv {output.tsv}
        """

def get_reference_assembly_targets(wildcards):
    f = checkpoints.select_references_meta.get().output.tsv
    df = pd.read_csv(f, sep="\\t")
    targets = []
    for i, row in df.iterrows():
        sample = row["sample"]
        segment = row["segment"]
        targets.append(config["output"] + f"assembly/{segment}/consensus/final_consensus/samples_alignment.fasta")
    return targets

rule run_reference_assemblies:
    input:
        get_reference_assembly_targets
    output:
        config["output"] + "assembly/reference_assemblies_done.txt"
    shell:
        "echo 'all dynamic reference assemblies complete' > {output}"

# Extract fasta for the given segment (family_accession) from the blast database
rule extract_reference_fasta:
    input:
        tsv = rules.select_references_meta.output.tsv
    output:
        fasta = config["output"] + "assembly/{segment}/references/{sample}.fasta"
    params:
        db = config.get("databases", {}).get("viral_genomes", "databases/virus_genomes/viral.genomes.fasta")
    conda:
        "../envs/utils.yaml"
    shell:
        """
        # Parse the TSV to get the accession for the segment/sample
        acc=$(awk -F'\\t' -v s="{wildcards.sample}" -v seg="{wildcards.segment}" 'NR>1 && $1==s && $2==seg {{print $3}}' {input.tsv} | head -n 1)
        if [ -n "$acc" ]; then
            # Extract sequence from fasta
            awk -v seq="$acc" 'BEGIN {{RS=">"; FS="\\n"}} $1~seq {{print ">"$0}}' {params.db} > {output.fasta}
        else
            touch {output.fasta}
        fi
        """

def get_meta_reference(wildcards):
    return config["output"] + f"assembly/{wildcards.segment}/references/{wildcards.sample}.fasta"

REFERENCE = get_meta_reference
SEGMENT_WILDCARD = "{segment}/"

if config["data_type"] == "illumina":
    include: "alignment_illumina.smk"
    include: "consensus_illumina.smk"
else:
    include: "alignment_nanopore.smk"
    include: "consensus_nanopore.smk"

