# Estrutura de Saída

## Pipeline de consenso (`viralunity consensus`)

O diretório de saída é organizado da seguinte forma após uma execução bem-sucedida:

```
{saida}/
├── assembly/
│   ├── assembly_stats_summary.csv    # métricas de QC por amostra
│   ├── consensus/                    # (não segmentado)
│   │   └── final_consensus/
│   │       └── samples_alignment.fasta
│   └── {segmento}/                  # (segmentado, um por segmento)
│       └── consensus/
│           └── final_consensus/
│               └── samples_alignment.fasta
├── reference/
│   └── reference.sanitized.fasta    # referência com cabeçalhos sanitizados (nanopore)
├── samples/
│   └── {nome_amostra}/
│       ├── consensus.fasta
│       ├── consensus.vcf.gz
│       ├── raw.vcf.gz
│       ├── table_cov_basewise.txt
│       ├── raw_mapped_reads.bam
│       └── trimmed_mapped_reads.bam
└── benchmark.tsv                    # tempo de execução por tarefa
```

### Arquivos principais

| Arquivo | Descrição |
|---------|-----------|
| `assembly/assembly_stats_summary.csv` | Profundidade de cobertura, cobertura horizontal, reads mapeadas por amostra/segmento |
| `samples/{amostra}/consensus.fasta` | Sequência de consenso final |
| `samples/{amostra}/table_cov_basewise.txt` | Tabela de cobertura por base |
| `samples/{amostra}/raw_mapped_reads.bam` | Reads mapeadas na referência |
| `benchmark.tsv` | Tempo de execução e uso de recursos por tarefa |

---

## Pipeline de metagenômica (`viralunity meta`)

```
{saida}/{run_name}/
├── qc/                                         # somente Illumina
│   ├── trim.{amostra}_fastp.html
│   ├── trim.{amostra}_fastp.json
│   └── reports/multiqc_report.html
├── host_filtered/                              # quando remoção do hospedeiro está ativa
│   ├── {amostra}.R1.filtered.fastq.gz
│   ├── {amostra}.R2.filtered.fastq.gz
│   └── {amostra}.merged.fastq.gz
├── metagenomics/
│   └── taxonomic_assignments/
│       ├── kraken2_reads/                      # quando --run-kraken2-reads
│       │   ├── results/{amostra}.report.txt
│       │   ├── results/{amostra}.output.krona.html
│       │   ├── kraken2_reads_taxa_summary.tsv
│       │   ├── kraken2_reads_taxa_summary_RPM.tsv
│       │   └── kraken2_reads_taxa_summary_RPM.bleed.tsv
│       ├── kraken2_contigs/                    # quando --run-denovo-assembly + --run-kraken2-contigs
│       │   ├── results/{amostra}.report.txt
│       │   ├── results/{amostra}.output.krona.html
│       │   └── kraken2_contigs_taxa_summary.tsv
│       ├── diamond_reads/                      # quando --run-diamond-reads
│       │   ├── results/{amostra}.diamond.tsv
│       │   ├── results/{amostra}.diamond.krona.html
│       │   ├── diamond_reads_taxa_summary.tsv
│       │   ├── diamond_reads_taxa_summary_RPM.tsv
│       │   └── diamond_reads_taxa_summary_RPM.bleed.tsv
│       └── diamond_contigs/                    # quando --run-denovo-assembly + --run-diamond-contigs
│           ├── results/{amostra}.diamond.supported.tsv
│           ├── results/{amostra}.diamond.supported.krona.html
│           ├── diamond_contigs_taxa_summary.tsv
│           └── diamond_contigs_taxa_summary_RPM.bleed.tsv
├── denovo_assembly/                            # quando --run-denovo-assembly
│   ├── megahit/{amostra}/final.contigs.fa
│   └── viral_contigs/{amostra}.viral_contigs.fa
├── assembly/                                   # quando --run-reference-assembly
│   └── {ref_key}/                              # chave única: {família}_{accession}
│       ├── references/{amostra}.fasta          # genoma de referência extraído
│       ├── mapped_reads/
│       │   ├── raw/{amostra}.sorted.bam
│       │   └── trimmed/{amostra}.sorted.bam
│       ├── isnvs/{amostra}.isnvs.vcf.gz        # somente Illumina, quando --run-isnv
│       └── consensus/final_consensus/
│           ├── {amostra}.consensus.fasta
│           └── {amostra}.consensus.vcf.gz
├── reference_targets.tsv                       # checkpoint: referências selecionadas por amostra/ref_key
├── reference_assembly_done.txt                 # sentinela: montagem de referência concluída
├── samples/                                    # links simbólicos por amostra para conveniência
│   └── {amostra}/
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
└── benchmark.tsv                               # tempo de execução e recursos por tarefa
```

### Arquivos principais

| Arquivo | Descrição |
|---------|-----------|
| `metagenomics/taxonomic_assignments/kraken2_reads/kraken2_reads_taxa_summary_RPM.bleed.tsv` | Tabela de táxons Kraken2 (reads) com normalização RPM e filtro de bleed |
| `metagenomics/taxonomic_assignments/diamond_reads/diamond_reads_taxa_summary_RPM.bleed.tsv` | Tabela de táxons Diamond (reads) com normalização RPM e filtro de bleed |
| `reference_targets.tsv` | Mapeia cada amostra × ref_key para o accession de referência selecionado |
| `assembly/{ref_key}/consensus/final_consensus/{amostra}.consensus.fasta` | Sequência de consenso guiada por referência por amostra e ref_key |
| `samples/{amostra}/` | Links simbólicos para todos os arquivos de saída por amostra |
| `benchmark.tsv` | Tempo de execução e uso de recursos por tarefa |
