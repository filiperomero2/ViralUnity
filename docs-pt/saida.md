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
├── qc/                              # somente Illumina
│   ├── {amostra}.fastp.json
│   ├── {amostra}.fastp.html
│   └── multiqc_report.html
├── host_filtered/                   # quando remoção do hospedeiro está ativa
│   └── {amostra}.host_filtered.fastq.gz
├── metagenomics/
│   └── taxonomic_assignments/
│       ├── {amostra}.kraken2.report
│       ├── {amostra}.krona.html
│       └── taxa_summary.tsv
├── denovo_assembly/                 # quando --run-denovo-assembly
│   └── {amostra}/
│       └── final.contigs.fa
└── medaka_work/                     # somente Nanopore, quando --run-polish-medaka
    └── {amostra}/
```
