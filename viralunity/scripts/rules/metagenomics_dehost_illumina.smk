if dehost_with_deacon:
    rule remove_host_reads:
        conda:
            "../envs/alignment.yaml"
        input:
            paired_R1 = rules.perform_qc.output.paired_R1,
            paired_R2 = rules.perform_qc.output.paired_R2,
            index = config["deacon_index"]
        output:
            filtered_R1 = config["output"] + "host_filtered/{sample}.R1.filtered.fastq.gz",
            filtered_R2 = config["output"] + "host_filtered/{sample}.R2.filtered.fastq.gz",
            summary = config["output"] + "logs/remove_host/{sample}.deacon_summary.json",
        threads: config["threads"]
        log:
            config["output"] + "logs/remove_host/{sample}.log"
        benchmark:
            config["output"] + "logs/remove_host/{sample}.benchmark.txt"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.filtered_R1}) $(dirname {log})
            deacon filter -d -t {threads} -s {output.summary} "{input.index}" "{input.paired_R1}" "{input.paired_R2}" \
                -o {output.filtered_R1} -O {output.filtered_R2} >> {log} 2>&1
            """

elif host_filtering_enabled:
    rule index_host_genome:
        conda:
            "../envs/utils.yaml"
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
            """
            set -euo pipefail
            minimap2 -d {output.index} {input.fasta} > {log} 2>&1
            """

    rule remove_host_reads:
        conda:
            "../envs/alignment.yaml"
        input:
            paired_R1 = rules.perform_qc.output.paired_R1,
            paired_R2 = rules.perform_qc.output.paired_R2,
            index = config["host_reference"] + ".mmi"
        output:
            filtered_R1 = config["output"] + "host_filtered/{sample}.R1.filtered.fastq.gz",
            filtered_R2 = config["output"] + "host_filtered/{sample}.R2.filtered.fastq.gz",
        threads: config["threads"]
        log:
            config["output"] + "logs/remove_host/{sample}.log"
        benchmark:
            config["output"] + "logs/remove_host/{sample}.benchmark.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.filtered_R1}) $(dirname {log})
            tmp_bam=$(mktemp --suffix .unmapped.bam)
            minimap2 -t {threads} -ax sr {input.index} {input.paired_R1} {input.paired_R2} \
            | samtools view -@ {threads} -b -f 4 - > "$tmp_bam"
            unmapped=$(samtools view -c "$tmp_bam")
            if [ "$unmapped" -eq 0 ]; then
                touch {output.filtered_R1} {output.filtered_R2}
                echo "All reads mapped to host; empty FASTQs written." >> {log}
            else
                samtools fastq -@ {threads} -1 {output.filtered_R1} -2 {output.filtered_R2} "$tmp_bam" >> {log} 2>&1
            fi
            rm -f "$tmp_bam"
            """

else:
    rule remove_host_reads:
        conda:
            "../envs/alignment.yaml"
        input:
            paired_R1 = rules.perform_qc.output.paired_R1,
            paired_R2 = rules.perform_qc.output.paired_R2,
        output:
            filtered_R1 = config["output"] + "host_filtered/{sample}.R1.filtered.fastq.gz",
            filtered_R2 = config["output"] + "host_filtered/{sample}.R2.filtered.fastq.gz",
        log:
            config["output"] + "logs/remove_host/{sample}.log"
        benchmark:
            config["output"] + "logs/remove_host/{sample}.benchmark.txt"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.filtered_R1}) $(dirname {log})
            cp {input.paired_R1} {output.filtered_R1}
            cp {input.paired_R2} {output.filtered_R2}
            echo "No host genome provided — skipping host read removal." > {log}
            """

rule merge_host_filtered_reads:
    conda:
        "../envs/utils.yaml"
    input:
        filtered_R1 = rules.remove_host_reads.output.filtered_R1,
        filtered_R2 = rules.remove_host_reads.output.filtered_R2,
    output:
        merged = temp(config["output"] + "host_filtered/{sample}.merged.fastq.gz"),
    log:
        config["output"] + "logs/merge_reads/{sample}.log"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.merged}) $(dirname {log})
        decompress() {{
            case "$1" in
                *.gz) gzip -dc "$1" 2>/dev/null || true ;;
                *) [ -s "$1" ] && cat "$1" || true ;;
            esac
        }}
        {{ decompress "{input.filtered_R1}"; decompress "{input.filtered_R2}"; }} | gzip -c > {output.merged}
        """
