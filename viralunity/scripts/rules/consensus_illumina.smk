# Illumina consensus rules: infer_consensus_sequence, detect_isnv, generate_vcf_consensus
# These rules expect the following variables to be defined in the entry-point workflow:
# - REFERENCE: path to reference genome (str)
# - config: standard Snakemake config dict


rule detect_isnv:
    conda:
        "../envs/alignment.yaml"
    input:
        reference = REFERENCE,
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        bam = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "isnvs/{sample}.lofreq.sorted.bam"),
        bam_index = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "isnvs/{sample}.lofreq.sorted.bam.bai"),
        vcf_tmp = temp(config['output'] + "assembly/" + SEGMENT_WILDCARD + "isnvs/{sample}.lofreq.tmp.vcf"),
        vcf = config['output'] + "assembly/" + SEGMENT_WILDCARD + "isnvs/{sample}.isnvs.vcf.gz",
        vcf_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "isnvs/{sample}.isnvs.vcf.gz.tbi"
    params:
        af_min_threshold = config["af_isnv_threshold"]
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/lofreq/{sample}.log"
    benchmark:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/lofreq/{sample}.benchmark.txt"
    threads: config["threads"]
    shell:
        """
        lofreq indelqual \
            -f {input.reference} \
            -o {output.bam} \
            --dindel \
            {input.bam}
        samtools index {output.bam} {output.bam_index}

        lofreq call-parallel \
            --pp-threads {threads} \
            --call-indels \
            -f {input.reference} \
            -o {output.vcf_tmp} \
            {output.bam}

        bcftools view -i 'INFO/AF<0.5 & INFO/AF>={params.af_min_threshold}' {output.vcf_tmp} -Oz -o {output.vcf}
        tabix {output.vcf}
        """

rule infer_consensus_sequence:
    conda:
        "../envs/alignment.yaml"
    input:
        bam = rules.trim_primer_sequences.output.bam,
        bam_index = rules.trim_primer_sequences.output.bam_index
    output:
        consensus = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.consensus.fasta"
    params:
        minimum_depth = config["minimum_depth"],
        af_threshold = config["af_threshold"]
    benchmark:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/samtools/consensus/{sample}.benchmark.txt"
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/samtools/consensus/{sample}.log"
    shell:
        "samtools consensus -a -d {params.minimum_depth} -m simple -q -c {params.af_threshold} --show-ins yes {input.bam} -o {output.consensus}"

rule generate_vcf_consensus:
    conda:
        "../envs/alignment.yaml"
    input:
        reference = REFERENCE,
        consensus = rules.infer_consensus_sequence.output.consensus
    output:
        vcf = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.consensus.vcf.gz",
        vcf_index = config['output'] + "assembly/" + SEGMENT_WILDCARD + "consensus/final_consensus/{sample}.consensus.vcf.gz.tbi"
    benchmark:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/gsaalign/{sample}.benchmark.txt"
    log:
        config['output'] + "assembly/" + SEGMENT_WILDCARD + "logs/gsaalign/{sample}.log"
    shell:
        """
        out_prefix=$(echo {output.vcf} | sed 's/.vcf.gz//')
        
        # Check if the consensus FASTA has sequence content (excluding header)
        if grep -v "^>" {input.consensus} | grep -q "[A-Za-z]"; then
            GSAlign \
                -r {input.reference} \
                -q {input.consensus} \
                -o $out_prefix \
                -fmt 1 \
                -sen

            bgzip -f $out_prefix.vcf
            tabix -p vcf {output.vcf} || touch {output.vcf_index}
            rm -f $out_prefix.maf
        else
            echo "Warning: Consensus sequence for {wildcards.sample} is empty. Creating a mock VCF." >&2
            echo "##fileformat=VCFv4.2" > $out_prefix.vcf
            printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" >> $out_prefix.vcf
            bgzip -f $out_prefix.vcf
            tabix -p vcf {output.vcf} || touch {output.vcf_index}
        fi
        """
