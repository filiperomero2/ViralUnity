if dehost_with_deacon:
    rule remove_host_reads:
        input:
            reads = lambda wildcards: config["samples"][wildcards.sample],
            index = config["deacon_index"]
        output:
            filtered = config["output"] + "host_filtered/{sample}.filtered.fastq.gz",
            summary = config["output"] + "logs/remove_host/{sample}.deacon_summary.json",
        threads: config.get("remove_host_reads_cpus", 2)
        resources:
            mem_mb = config.get("remove_host_reads_ram", 4) * 1024
        log:
            config["output"] + "logs/remove_host/{sample}.log"
        benchmark:
            config["output"] + "logs/remove_host/{sample}.benchmark.txt"
        conda:
            "../envs/alignment.yaml"
        shell:
            """
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
        threads: config.get("index_host_genome_cpus", 2)
        resources:
            mem_mb = config.get("index_host_genome_ram", 4) * 1024
        log:
            config["output"] + "logs/remove_host/host_genome_indexing.log"
        benchmark:
            config["output"] + "logs/remove_host/host_genome_indexing.benchmark.txt"
        conda:
            "../envs/utils.yaml"
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
        threads: config.get("remove_host_reads_cpus", 2)
        resources:
            mem_mb = config.get("remove_host_reads_ram", 4) * 1024
        log:
            config["output"] + "logs/remove_host/{sample}.log"
        benchmark:
            config["output"] + "logs/remove_host/{sample}.benchmark.txt"
        conda:
            "../envs/alignment.yaml"
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
        conda:
            "../envs/alignment.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {output.filtered})" "$(dirname {log})"
            if [[ "{input.reads}" == *.gz ]]; then
                cp -f "{input.reads}" "{output.filtered}"
            else
                gzip -c "{input.reads}" > "{output.filtered}"
            fi
            echo "No host genome provided — skipping host read removal." > "{log}"
            """
