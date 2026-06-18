rule merge_nanopore_reads:
    input:
        reads = get_map_input_fastqs
    output:
        merged = temp(config["output"] + "raw_merged/{sample}.fastq.gz")
    log:
        config["output"] + "logs/merge_reads/{sample}.log"
    benchmark:
        config["output"] + "logs/merge_reads/{sample}.benchmark.txt"
    conda:
        "../envs/utils.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.merged})" "$(dirname {log})"
        is_read_file() {{
            case "$(basename "$1")" in
                *[Zz][Oo][Nn][Ee].[Ii]dentifier*) return 1 ;;
                *.fastq|*.fastq.gz|*.fq|*.fq.gz|*.fasta|*.fasta.gz) return 0 ;;
                *) return 1 ;;
            esac
        }}
        reads=()
        for _f in {input.reads}; do
            if is_read_file "$_f"; then
                reads+=("$_f")
            fi
        done
        if [ "${{#reads[@]}}" -eq 0 ]; then
            echo "No read files left after filtering non-FASTQ paths (e.g. Zone.Identifier)." >&2
            exit 1
        fi
        if [ "${{#reads[@]}}" -eq 1 ]; then
            f="${{reads[0]}}"
            if [[ "$f" == *.gz ]]; then
                cp -f "$f" "{output.merged}"
            else
                gzip -c "$f" > "{output.merged}"
            fi
        else
            decompress() {{
                case "$1" in
                    *.gz) gzip -dc "$1" 2>/dev/null || true ;;
                    *) [ -s "$1" ] && cat "$1" || true ;;
                esac
            }}
            for f in "${{reads[@]}}"; do
                decompress "$f"
            done | gzip -c > "{output.merged}"
        fi
        """

if dehost_with_deacon:
    rule remove_host_reads:
        input:
            reads = rules.merge_nanopore_reads.output.merged,
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
            "../envs/alignment.yaml"
        shell:
            r"""
            set -euo pipefail
            minimap2 -d {output.index} {input.fasta} 2> {log}
            """

    rule remove_host_reads:
        input:
            reads = rules.merge_nanopore_reads.output.merged,
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
            reads = rules.merge_nanopore_reads.output.merged,
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
