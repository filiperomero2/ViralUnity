# Nanopore consensus rules: infer_consensus_sequence (clair3 + bcftools)
# These rules expect the following variables to be defined in the entry-point workflow:
# - REFERENCE: path to reference genome (str)
# - config: standard Snakemake config dict


rule infer_consensus_sequence:
    conda:
        "../envs/clair3.yaml"
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index,
        reference = REFERENCE
    output:
        vcf_raw = config['output'] + "assembly/" + SEGMENT_WILDCARD + "clair3/{sample}/{sample}.raw.vcf.gz",
        vcf_raw_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "clair3/{sample}/{sample}.raw.vcf.gz.tbi",
        vcf_norm = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.norm.vcf.gz"),
        vcf = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.vcf.gz",
        vcf_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.vcf.gz.tbi",
        low_cov_bed = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.low_cov.bed",
        consensus = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.consensus.fasta"
    params:
        output_prefix_dir = config['output'] + "assembly/" + SEGMENT_WILDCARD + "clair3/{sample}",
        minimum_depth = config.get("minimum_depth", 10),
        af_threshold = config.get("af_threshold", 0.7),
        chunk_size = config.get("chunk_size", 50000),
        clair3_model = config.get("clair3_model", "r1041_e82_400bps_sup_v420"),
        variant_quality = config.get("variant_quality", 20),
        minimum_map_quality = config.get("minimum_map_quality", 20),
        variant_depth = config.get("variant_depth", 5)
    benchmark:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/consensus/{sample}.benchmark.txt"
    threads: config.get("infer_consensus_sequence_cpus", 2)
    resources:
        mem_mb = config.get("infer_consensus_sequence_ram", 4) * 1024
    shell:
        """
        set -euo pipefail
        run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.reference} \
            --qual={params.variant_quality} \
            --min_mq={params.minimum_map_quality} \
            --model_path=$CONDA_PREFIX/bin/models/{params.clair3_model} \
            --chunk_size={params.chunk_size} \
            --threads={threads} \
            --enable_long_indel \
            --haploid_sensitive \
            --no_phasing_for_fa \
            --output={params.output_prefix_dir} \
            --platform='ont' \
            --include_all_ctgs 
        
        cp {params.output_prefix_dir}/merge_output.vcf.gz {output.vcf_raw}
        cp {params.output_prefix_dir}/merge_output.vcf.gz.tbi {output.vcf_raw_index}

        bcftools norm -m - -f {input.reference} {output.vcf_raw} > {output.vcf_norm}
        bcftools filter -i 'FILTER="PASS" && FORMAT/AF >= {params.af_threshold} && FORMAT/AD[0:1] >= {params.variant_depth}' {output.vcf_norm} -o {output.vcf} -O z
        tabix {output.vcf}
        
        samtools depth -J -a {input.bam} | \
            awk '$3 <= int({params.minimum_depth}) {{print $1 "\t" $2-1 "\t" $2}}' > {output.low_cov_bed}
        
        bcftools consensus -f {input.reference} --mask {output.low_cov_bed} {output.vcf} > {output.consensus}

        find {params.output_prefix_dir} -mindepth 1 ! -name '*.raw.vcf.gz' ! -name '*.raw.vcf.gz.tbi' -exec rm -rf {{}} +
        """
