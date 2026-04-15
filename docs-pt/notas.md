# Notas

## Vírus segmentados

O ViralUnity suporta nativamente a montagem de genomas virais segmentados. Em vez de executar o pipeline várias vezes, passe múltiplas referências de segmentos com `--segmented-reference`:

```bash
viralunity consensus illumina \
    --sample-sheet amostras.csv \
    --config-file config_segmentado.yml \
    --output /caminho/saida \
    --segmented-reference S=/caminho/segmento_s.fasta \
    --segmented-reference L=/caminho/segmento_l.fasta
```

Quando esse argumento é detectado, o pipeline aciona um workflow modular especializado que processa cada segmento independentemente em paralelo. Os resultados são organizados em `samples/{nome_amostra}/{nome_segmento}/`.

## Dados Nanopore

O pipeline de metagenômica suporta dados Nanopore via `viralunity meta nanopore`. Não há etapa de QC com fastp; remoção opcional do hospedeiro, montagem com MEGAHIT e polimento com Racon/Medaka estão disponíveis. As planilhas de amostras possuem duas colunas (ID da amostra, caminho para um FASTQ/FASTA por amostra).

## Sanitização de cabeçalhos de referência

O workflow nanopore de consenso sanitiza automaticamente os cabeçalhos do FASTA de referência antes do uso. Caracteres especiais (`/`, `\`, `|`, `,`, `~` e espaços) nos identificadores de sequências são substituídos por underscores (`_`). Isso previne problemas em ferramentas posteriores, como o clair3, que utiliza o ID da sequência para criar diretórios de saída.

## Montagem dinâmica por referência

Quando `--run-reference-assembly` está habilitado, o pipeline de metagenômica adiciona uma etapa pós-classificação que seleciona automaticamente genomas de referência e executa montagem de consenso para cada amostra que atinge o limiar de hits.

### Como funciona

1. **Checkpoint de seleção** — após a classificação taxonômica, um checkpoint do Snakemake lê os arquivos de resumo de táxons e seleciona os candidatos a referência. As flags `--method` e `--source` controlam quais resultados de classificação são utilizados (ex: `kraken2` em `reads`). Apenas amostras com pelo menos `--reads-count` reads (ou `--contigs-count` contigs) atribuídas a uma das `--families` são selecionadas.

2. **Extração de referência** — para cada amostra selecionada, uma chave de montagem única (`ref_key`) é construída a partir da família viral e do accession de referência (`{família}_{accession}`). O genoma correspondente é extraído do arquivo `--viral-genomes`. Duas estratégias de seleção estão disponíveis:
   - `taxid` (padrão): o taxid do resultado de classificação é consultado em `--viral-taxids` para encontrar o accession correspondente.
   - `similarity`: os contigs de novo são alinhados via BLAST contra `--viral-genomes` e o melhor hit acima dos limiares `--blast-qcov` / `--blast-pident` é utilizado.

3. **Montagem de consenso** — utiliza as mesmas regras de alinhamento e consenso do `viralunity consensus`, com a referência selecionada dinamicamente. Os resultados são gerados em `assembly/{ref_key}/`.

### Bancos de dados necessários

| Opção | Arquivo gerado por |
|-------|--------------------|
| `--viral-genomes` | `viralunity get-databases virus-genome` |
| `--viral-taxids` | `viralunity get-databases virus-genome` |

### Exemplo

```bash
viralunity meta illumina \
    --sample-sheet amostras.csv \
    --config-file config.yaml \
    --output /caminho/saida \
    --kraken2-database /caminho/kraken2_db \
    --krona-database /caminho/krona_taxonomy \
    --taxdump /caminho/taxdump \
    --run-reference-assembly \
    --method kraken2 \
    --source reads \
    --families Coronaviridae,Orthomyxoviridae \
    --reads-count 100 \
    --viral-genomes databases/virus_genomes/viral.genomes.fasta \
    --viral-taxids databases/virus_genomes/genome2taxid.tsv \
    --threads 4 --threads-total 8
```

## Executar os testes

```bash
pytest -v test
```
