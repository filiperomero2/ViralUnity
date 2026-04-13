if run_denovo:

    rule run_megahit:
        input:
            reads = rules.remove_host_reads.output.filtered
        output:
            contigs = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa"
        threads: config.get("run_megahit_cpus", 2)
        resources:
            mem_mb = config.get("run_megahit_ram", 4) * 1024
        log:
            config["output"] + "logs/megahit/{sample}.log"
        benchmark:
            config["output"] + "logs/megahit/{sample}.benchmark.txt"
        params:
            tempdir = config["output"] + "denovo_assembly/megahit/temp_{sample}",
            outdir = config["output"] + "denovo_assembly/megahit/{sample}"
        conda:
            "../envs/assembly.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {log})
            workdir="${{TMPDIR:-/tmp}}/megahit_{wildcards.sample}"
            mkdir -p "$workdir"
            if [ ! -s {input.reads} ]; then
                echo "Reads empty; skipping MEGAHIT." > {log}
                mkdir -p $(dirname {output.contigs})
                touch {output.contigs}
            else
                megahit -r {input.reads} -o {params.tempdir} --num-cpu-threads {threads} --tmp-dir "$workdir" >> {log} 2>&1
                mv {params.tempdir}/final.contigs.fa {output.contigs}
                rm -rf {params.tempdir}
            fi
            """

    if run_polish_racon:

        rule run_racon:
            input:
                assembly = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa",
                reads = rules.remove_host_reads.output.filtered
            output:
                racon_fasta = config["output"] + "denovo_assembly/megahit/{sample}/racon.fasta"
            threads: config.get("run_racon_cpus", 2)
            resources:
                mem_mb = config.get("run_racon_ram", 4) * 1024
            log:
                config["output"] + "logs/racon/{sample}.log"
            benchmark:
                config["output"] + "logs/racon/{sample}.benchmark.txt"
            conda:
                "../envs/consensus.yaml"
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {log})
                if [ ! -s {input.assembly} ]; then
                    echo "Assembly empty; skipping Racon." > {log}
                    touch {output.racon_fasta}
                else
                    workdir="${{TMPDIR:-/tmp}}/racon_{wildcards.sample}"
                    mkdir -p "$workdir"
                    minimap2 -t {threads} -x map-ont {input.assembly} {input.reads} > "$workdir/align.paf" 2>> {log}
                    racon -t {threads} {input.reads} "$workdir/align.paf" {input.assembly} > {output.racon_fasta} 2>> {log}
                fi
                """

    if run_polish_medaka:

        rule run_medaka:
            input:
                assembly = get_medaka_assembly_input,
                reads = rules.remove_host_reads.output.filtered
            output:
                polished = config["output"] + "denovo_assembly/megahit/{sample}/polished.fasta",
                bam = config["output"] + "medaka_work/{sample}/calls_to_draft.bam"
            threads: config.get("run_medaka_cpus", 2)
            resources:
                mem_mb = config.get("run_medaka_ram", 4) * 1024
            log:
                config["output"] + "logs/medaka/{sample}.log"
            benchmark:
                config["output"] + "logs/medaka/{sample}.benchmark.txt"
            params:
                outdir = config["output"] + "medaka_work/{sample}"
            conda:
                "../envs/medaka.yaml"
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {log})
                if [ ! -s {input.assembly} ]; then
                    echo "No assembly to polish for sample {wildcards.sample}." > {log}
                    mkdir -p $(dirname {output.polished})
                    touch {output.polished}
                    mkdir -p {params.outdir}
                    touch {output.bam}
                else
                    mkdir -p {params.outdir}
                    medaka_consensus \
                        -i {input.reads} \
                        -d {input.assembly} \
                        -o {params.outdir} \
                        -g -r 'N' \
                        -t {threads} &> {log}
                    mv {params.outdir}/consensus.fasta {output.polished}
                fi
                """
