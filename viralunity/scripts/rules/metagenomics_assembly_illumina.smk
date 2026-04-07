if run_denovo:

    rule run_megahit:
        input:
            filtered_R1 = rules.remove_host_reads.output.filtered_R1,
            filtered_R2 = rules.remove_host_reads.output.filtered_R2,
        output:
            contigs = config["output"] + "denovo_assembly/megahit/{sample}/final.contigs.fa"
        threads: config["threads"]
        log:
            config["output"] + "logs/megahit/{sample}.log"
        benchmark:
            config["output"] + "logs/megahit/{sample}.benchmark.txt"
        params:
            tempdir = temp(config["output"] + "denovo_assembly/megahit/temp_{sample}"),
            outdir = config["output"] + "denovo_assembly/megahit/{sample}"
        conda:
            "../envs/assembly.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {log})
            workdir="${{TMPDIR:-/tmp}}/megahit_{wildcards.sample}"
            mkdir -p "$workdir"
            size_R1=$([ -s {input.filtered_R1} ] && echo 1 || echo 0)
            size_R2=$([ -s {input.filtered_R2} ] && echo 1 || echo 0)
            if [ "$size_R1" -eq 1 ] && [ "$size_R2" -eq 1 ]; then
                megahit -1 {input.filtered_R1} -2 {input.filtered_R2} -o {params.tempdir} \
                    --num-cpu-threads {threads} --tmp-dir "$workdir" >> {log} 2>&1
                mv {params.tempdir} {params.outdir}
            else
                echo "R1 and/or R2 empty; skipping MEGAHIT." > {log}
                mkdir -p $(dirname {output.contigs})
                touch {output.contigs}
            fi
            """
