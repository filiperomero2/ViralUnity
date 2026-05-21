# Output Layout

## Consensus pipeline (`viralunity consensus`)

After a successful run, the output directory is organised as follows:

```
{output}/
├── assembly/
│   ├── assembly_stats_summary.csv    # per-sample QC metrics
│   ├── consensus/                    # non-segmented
│   │   └── final_consensus/
│   │       └── samples_alignment.fasta
│   └── {segment}/                    # segmented (one per segment)
│       └── consensus/
│           └── final_consensus/
│               └── samples_alignment.fasta
├── reference/
│   └── reference.sanitized.fasta    # reference with sanitised headers (nanopore)
├── samples/
│   └── {sample_name}/
│       ├── consensus.fasta
│       ├── consensus.vcf.gz
│       ├── raw.vcf.gz
│       ├── table_cov_basewise.txt
│       ├── raw_mapped_reads.bam
│       └── trimmed_mapped_reads.bam
└── benchmark.tsv                    # per-task runtime
```

### Key files

| File | Description |
|------|-------------|
| `assembly/assembly_stats_summary.csv` | Depth of coverage, breadth of coverage, mapped reads per sample/segment |
| `samples/{sample}/consensus.fasta` | Final consensus sequence |
| `samples/{sample}/table_cov_basewise.txt` | Per-base coverage table |
| `samples/{sample}/raw_mapped_reads.bam` | Reads mapped to the reference |
| `benchmark.tsv` | Runtime and resource usage per task |

---

## Metagenomics pipeline (`viralunity meta`)

```
{output}/{run_name}/
├── qc/                                         # Illumina only
│   ├── trim.{sample}_fastp.html
│   ├── trim.{sample}_fastp.json
│   └── reports/multiqc_report.html
├── host_filtered/                              # when host removal is enabled
│   ├── {sample}.R1.filtered.fastq.gz
│   ├── {sample}.R2.filtered.fastq.gz
│   └── {sample}.merged.fastq.gz
├── metagenomics/
│   └── taxonomic_assignments/
│       ├── kraken2_reads/                      # when --run-kraken2-reads
│       │   ├── results/{sample}.report.txt
│       │   ├── results/{sample}.output.krona.html
│       │   ├── kraken2_reads_taxa_summary.tsv
│       │   ├── kraken2_reads_taxa_summary_RPM.tsv
│       │   └── kraken2_reads_taxa_summary_RPM.bleed.tsv
│       ├── kraken2_contigs/                    # when --run-denovo-assembly + --run-kraken2-contigs
│       │   ├── results/{sample}.report.txt
│       │   ├── results/{sample}.output.krona.html
│       │   └── kraken2_contigs_taxa_summary.tsv
│       ├── diamond_reads/                      # when --run-diamond-reads
│       │   ├── results/{sample}.diamond.tsv
│       │   ├── results/{sample}.diamond.krona.html
│       │   ├── diamond_reads_taxa_summary.tsv
│       │   ├── diamond_reads_taxa_summary_RPM.tsv
│       │   └── diamond_reads_taxa_summary_RPM.bleed.tsv
│       └── diamond_contigs/                    # when --run-denovo-assembly + --run-diamond-contigs
│           ├── results/{sample}.diamond.supported.tsv
│           ├── results/{sample}.diamond.supported.krona.html
│           ├── diamond_contigs_taxa_summary.tsv
│           └── diamond_contigs_taxa_summary_RPM.bleed.tsv
├── denovo_assembly/                            # when --run-denovo-assembly
│   ├── megahit/{sample}/final.contigs.fa
│   └── viral_contigs/{sample}.viral_contigs.fa
├── assembly/                                   # when --run-reference-assembly
│   └── {ref_key}/                              # unique key: {family}_{accession}
│       ├── references/{sample}.fasta           # extracted reference genome
│       ├── mapped_reads/
│       │   ├── raw/{sample}.sorted.bam
│       │   └── trimmed/{sample}.sorted.bam
│       ├── isnvs/{sample}.isnvs.vcf.gz         # Illumina only, when --run-isnv
│       └── consensus/final_consensus/
│           ├── {sample}.consensus.fasta
│           └── {sample}.consensus.vcf.gz
├── reference_targets.tsv                       # checkpoint: selected references per sample/ref_key
├── reference_assembly_done.txt                 # sentinel: reference assembly completed
├── samples/                                    # per-sample symlinks for convenience
│   └── {sample}/
│       ├── fastp.html
│       ├── host_filtered_R1.fastq.gz
│       ├── host_filtered_R2.fastq.gz
│       ├── kraken2_reads.report.txt
│       ├── kraken2_reads.krona.html
│       ├── diamond_reads.tsv
│       ├── denovo_contigs.fasta
│       ├── kraken2_contigs.report.txt
│       ├── diamond_contigs_supported.tsv
│       └── viral_mapped_reads.bam
└── benchmark.tsv                               # runtime and resources per task
```

### Key files

| File | Description |
|------|-------------|
| `metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.tsv` | Kraken2 (reads) taxa table with RPM normalisation and bleed filter |
| `metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.tsv` | DIAMOND (reads) taxa table with RPM normalisation and bleed filter |
| `reference_targets.tsv` | Maps each sample × ref_key to the selected reference accession |
| `assembly/{ref_key}/consensus/final_consensus/{sample}.consensus.fasta` | Reference-guided consensus sequence per sample and ref_key |
| `samples/{sample}/` | Symlinks to every per-sample output for convenience |
| `benchmark.tsv` | Runtime and resource usage per task |
