# Experimental

rule all:
    input:
        config['output'] + "metagenomics/metagenomics_summary.txt"

def get_map_input_fastqs(wildcards):
    reads = config["samples"][wildcards.sample].split()
    return(reads)

def get_classification_filter():
    if config['remove_human_reads'] and config['remove_unclassified_reads']:
        my_filter = "\t9606|\t0"
    if not config['remove_human_reads'] and config['remove_unclassified_reads']:
        my_filter = "\t0"
    if config['remove_human_reads'] and not config['remove_unclassified_reads']:
        my_filter = "\t9606"
    if not config['remove_human_reads'] and not config['remove_unclassified_reads']:
        my_filter = " "
    print(my_filter)
    return(my_filter)

my_filter = get_classification_filter()     

rule run_kraken2:
    input:
        get_map_input_fastqs
    output:
        report = config['output'] + "metagenomics/taxonomic_assignments/results/{sample}.report.txt",
        outfile = config['output'] + "metagenomics/taxonomic_assignments/results/{sample}.output.txt",
    threads: 2
    params:
        database = config['kraken2_database']
    log:
        config['output'] + "logs/kraken2/{sample}.log"
    benchmark:
        config['output'] + "logs/kraken2/{sample}.benchmark.txt"
    shell:
        "kraken2 --db {params.database} --threads {threads} --report-minimizer-data "
        "--minimum-hit-group 3 --report {output.report} "
        "--output {output.outfile} {input} 2> {log}"

rule clean_kraken2_output:
    input:
        config['output'] + "metagenomics/taxonomic_assignments/results/{sample}.output.txt"
    output:
        config['output'] + "metagenomics/taxonomic_assignments/results/{sample}.output.krona.txt"
    params:
        to_rm = my_filter
    shell:
        "cat {input} | cut -f 2,3 | egrep -v \"{params.to_rm}\" > {output} "

rule create_krona_report:
    input:
        config['output'] + "metagenomics/taxonomic_assignments/results/{sample}.output.krona.txt"
    output:
        config['output'] + "metagenomics/taxonomic_assignments/reports/{sample}.output.krona.html"
    params:
        krona_database = config['krona_database']
    log:
        config['output'] + "logs/krona/{sample}.log"
    benchmark:
        config['output'] + "logs/krona/{sample}.benchmark.txt"
    shell:
        "ktImportTaxonomy {input} -tax {params.krona_database} -o {output} 2> {log}"

rule summarize:
    input:
        expand(config['output'] + "metagenomics/taxonomic_assignments/reports/{sample}.output.krona.html", sample=config["samples"])
    output:
        config['output'] + "metagenomics/metagenomics_summary.txt"
    params:
        path = config['output'] + "metagenomics/taxonomic_assignments/results/"
    shell:
        "echo \"file,percentage_of_reads,number_of_reads_rooted,number_of_reads_direct,number_of_k-mers,number_of_distinct_k-mers,rank_code,ncbi_taxid,taxon\" > {output}; "
        "grep \"viridae\" {params.path}*report.txt | sed -E \"s/[	| |:]+/,/g\">> {output}"