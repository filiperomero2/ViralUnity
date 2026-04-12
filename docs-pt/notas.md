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

## Executar os testes

```bash
pytest -v test
```
