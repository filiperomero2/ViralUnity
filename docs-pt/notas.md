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

2. **Extração de referência** — para cada amostra selecionada, uma chave de montagem única (`ref_key`) é construída a partir da família viral e do accession de referência (`{família}_{accession}`). O genoma correspondente é extraído do arquivo `--viral-genomes`. A estratégia de seleção controla como o accession de referência é escolhido (ver abaixo).

3. **Montagem de consenso** — utiliza as mesmas regras de alinhamento e consenso do `viralunity consensus`, com a referência selecionada dinamicamente. Os resultados são gerados em `assembly/{ref_key}/`.

### Estratégias de seleção de referência

#### `--reference-selection-strategy taxid` (padrão)

O taxid atribuído a cada read ou contig pelo classificador (Kraken2 ou Diamond) é consultado no arquivo `--viral-taxids` (`genome2taxid.tsv`) para encontrar todos os accessions de genoma que compartilham aquele taxid. Cada accession encontrado gera uma montagem de referência separada para aquela amostra.

**Quando usar:** quando o banco de dados do classificador e o banco `--viral-genomes` foram construídos a partir do mesmo release do RefSeq. A ligação pelo taxid é direta e não requer comparação de sequências. Rápida e determinística.

**Limitação:** uma correspondência de taxid não garante que a referência seja a mais próxima geneticamente — apenas garante identidade taxonômica. Para famílias muito diversas (ex.: Flaviviridae), múltiplos accessions podem corresponder ao mesmo taxid, gerando uma montagem para cada um.

**Bancos de dados necessários:** `--viral-genomes` e `--viral-taxids` (ambos gerados por `viralunity get-databases virus-genome`).

```bash
viralunity meta illumina \
    --sample-sheet amostras.csv \
    --config-file config.yaml \
    --output /caminho/saida \
    --kraken2-database /caminho/kraken2_db \
    --krona-database /caminho/krona_taxonomy \
    --taxdump /caminho/taxdump \
    --run-reference-assembly \
    --reference-selection-strategy taxid \
    --method kraken2 \
    --source reads \
    --families Coronaviridae,Orthomyxoviridae \
    --reads-count 100 \
    --viral-genomes databases/virus_genomes/viral.genomes.fasta \
    --viral-taxids databases/virus_genomes/genome2taxid.tsv \
    --threads 4 --threads-total 8
```

#### `--reference-selection-strategy similarity`

Os contigs montados de novo de cada amostra são alinhados via `blastn` contra `--viral-genomes`. O melhor hit que superar os limiares `--blast-qcov` e `--blast-pident` é selecionado como referência.

**Quando usar:** quando a identidade de sequência com um genoma específico é mais importante que o rótulo taxonômico — por exemplo, ao trabalhar com cepas divergentes ou novas, onde a busca por taxid pode retornar uma referência geneticamente distante. Requer montagem de novo (`--run-denovo-assembly`) para gerar os contigs usados como query no BLAST.

**Limitação:** requer um índice BLAST pré-construído ao lado do arquivo `--viral-genomes` (criado automaticamente por `viralunity get-databases virus-genome`). Se nenhum contig atingir os limiares de identidade/cobertura, nenhuma montagem de referência é iniciada para aquela amostra.

**Bancos de dados necessários:** `--viral-genomes` com seu índice BLAST (gerado por `viralunity get-databases virus-genome`).

```bash
viralunity meta illumina \
    --sample-sheet amostras.csv \
    --config-file config.yaml \
    --output /caminho/saida \
    --kraken2-database /caminho/kraken2_db \
    --krona-database /caminho/krona_taxonomy \
    --taxdump /caminho/taxdump \
    --run-denovo-assembly \
    --run-reference-assembly \
    --reference-selection-strategy similarity \
    --method kraken2 \
    --source reads \
    --families Coronaviridae,Orthomyxoviridae \
    --reads-count 100 \
    --blast-qcov 80 \
    --blast-pident 80 \
    --viral-genomes databases/virus_genomes/viral.genomes.fasta \
    --threads 4 --threads-total 8
```

### Comparação entre estratégias

| | `taxid` | `similarity` |
|---|---|---|
| Requer montagem de novo | Não | Sim |
| Requer `--viral-taxids` | Sim | Não |
| Requer índice BLAST | Não | Sim |
| Base da seleção | Correspondência de taxid | Identidade de sequência |
| Velocidade | Rápida | Mais lenta (BLAST por amostra) |
| Melhor para | Cepas conhecidas com boa cobertura no RefSeq | Cepas divergentes ou novas |

### Bancos de dados necessários

| Opção | Arquivo gerado por | Usado por |
|-------|--------------------|-----------| 
| `--viral-genomes` | `viralunity get-databases virus-genome` | Ambas as estratégias |
| `--viral-taxids` | `viralunity get-databases virus-genome` | Somente `taxid` |
| Índice BLAST (`.nhr`/`.nin`/`.nsq`) | `viralunity get-databases virus-genome` | Somente `similarity` |

## Executar os testes

```bash
pytest -v test
```
