# Referência de Comandos

## `viralunity create-samplesheet`

Antes de executar qualquer pipeline, você precisa de uma planilha CSV de amostras. O comando `create-samplesheet` a gera automaticamente a partir de um diretório de run.

```bash
viralunity create-samplesheet --input /caminho/para/run/ --output amostras.csv
```

### Opções

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--input` | *(obrigatório)* | Diretório do run: arquivos FASTQ diretamente (nível 0) ou um subdiretório por amostra (nível 1). |
| `--output` | *(obrigatório)* | Caminho do arquivo CSV de saída. |
| `--level` | `1` | `0` = arquivos em `--input`, `1` = arquivos em subdiretórios (um por amostra). |
| `--separator` | `-` | Caractere usado para dividir o nome do arquivo/diretório e extrair o ID da amostra (`-`, `_` ou `.`). |
| `--pattern` | `R1` | Padrão identificando o primeiro arquivo de leitura no nível 0 (`R1` para Illumina, `barcode` para Nanopore). |

### Formato de saída

O CSV gerado não possui linha de cabeçalho:

```
# Illumina (3 colunas)
amostra1,/caminho/amostra1_R1.fastq.gz,/caminho/amostra1_R2.fastq.gz

# Nanopore (2 colunas)
barcode01,/caminho/barcode01.fastq
```

### Exemplos

**Illumina** — FASTQs organizados em um subdiretório por amostra (nível 1, padrão):

```bash
viralunity create-samplesheet \
    --input /caminho/illumina_run/ \
    --output /caminho/amostras.csv
```

**Illumina** — todos os FASTQs diretamente no diretório do run (nível 0):

```bash
viralunity create-samplesheet \
    --input /caminho/illumina_run/ \
    --output /caminho/amostras.csv \
    --level 0 \
    --separator _ \
    --pattern R1
```

**Nanopore** — um FASTQ por subdiretório de barcode:

```bash
viralunity create-samplesheet \
    --input /caminho/nanopore_run/ \
    --output /caminho/amostras_nano.csv \
    --separator _ \
    --pattern barcode
```

---

## `viralunity consensus`

O pipeline de consenso leva reads brutas até sequências de consenso processadas com um único comando. Selecione o tipo de dados como subcomando (`illumina` ou `nanopore`).

### Opções — compartilhadas (ambos os tipos de dados)

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--sample-sheet` | *(obrigatório)* | Arquivo CSV com IDs de amostras e caminhos de arquivos. |
| `--config-file` | *(obrigatório)* | Caminho para o arquivo de configuração YAML a ser criado. |
| `--output` | *(obrigatório)* | Diretório base de saída. |
| `--reference` | — | FASTA de referência (mutuamente exclusivo com `--segmented-reference`). |
| `--segmented-reference` | — | Referência por segmento: `SEGMENTO=CAMINHO` (repetível). |
| `--primer-scheme` | — | Arquivo BED do esquema de primers (somente amplicon). |
| `--minimum-coverage` | `20` | Profundidade mínima para inclusão de base no consenso. |
| `--minimum-read-length` | `50` | Comprimento mínimo de read. |
| `--af-threshold` | `0.51` | Frequência alélica mínima para incluir variante no consenso. |
| `--run-name` | `undefined` | Nome do run de sequenciamento. |
| `--threads` | `1` | Threads por tarefa individual. |
| `--threads-total` | `1` | Total de threads para o workflow. |
| `--create-config-only` | off | Apenas gerar o arquivo de configuração; não executar o workflow. |

### Opções — somente Illumina

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--adapters` | — | FASTA de sequências adaptadoras (fastp QC). |
| `--trim-head` | `0` | Bases para remover do extremo 5′. |
| `--trim-tail` | `0` | Bases para remover do extremo 3′. |
| `--cut-front-mean-quality` | `10` | Limiar de qualidade cut_front do fastp. |
| `--cut-tail-mean-quality` | `10` | Limiar de qualidade cut_tail do fastp. |
| `--cut-right-window-size` | `4` | Tamanho da janela cut_right do fastp. |
| `--cut-right-mean-quality` | `15` | Limiar de qualidade cut_right do fastp. |
| `--af-isnv-threshold` | `0.0` | Frequência alélica mínima para análise de iSNV. |
| `--run-isnv` | off | Executar análise de SNV intra-hospedeiro (LoFreq). |

### Opções — somente Nanopore

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--chunk-size` | `10000` | Tamanho de chunk para processamento com clair3. |
| `--clair3-model` | `r1041_e82_400bps_sup_v500` | Modelo clair3 para chamada de variantes. |
| `--variant-quality` | `20` | Qualidade mínima de variante (clair3). |
| `--variant-depth` | `10` | Profundidade mínima do alelo alternativo (clair3). |
| `--minimum-map-quality` | `30` | Qualidade mínima de mapeamento (clair3). |

### Exemplos

**Illumina — referência única:**

```bash
viralunity consensus illumina \
    --sample-sheet /caminho/amostras.csv \
    --config-file /caminho/config.yml \
    --run-name exemplo \
    --output /caminho/saida \
    --reference /caminho/referencia.fasta \
    --threads 2 \
    --threads-total 4
```

**Illumina — genoma segmentado:**

```bash
viralunity consensus illumina \
    --sample-sheet /caminho/amostras.csv \
    --config-file /caminho/config_seg.yml \
    --run-name exemplo \
    --output /caminho/saida \
    --segmented-reference L=/caminho/segmento_L.fasta \
    --segmented-reference S=/caminho/segmento_S.fasta \
    --threads 2 \
    --threads-total 4
```

**Nanopore:**

```bash
viralunity consensus nanopore \
    --sample-sheet /caminho/amostras_nano.csv \
    --config-file /caminho/config_nano.yml \
    --run-name exemplo \
    --output /caminho/saida \
    --reference /caminho/referencia.fasta \
    --threads 4 \
    --threads-total 4
```

---

## `viralunity meta`

O pipeline de metagenômica tem como objetivo gerar classificações taxonômicas e visualizações com krona. Você pode usar **somente Kraken2**, **somente Diamond** ou **ambos**. Cada ferramenta pode ser executada em **reads** e, quando a montagem está habilitada, em **contigs**. No final da análise, é possível realizar montagem de consenso baseada em genomas de referência baixados sob demanda com base nos resultados das classificações ("dynamic reference assembly").

### Opções — compartilhadas (ambos os tipos de dados)

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--sample-sheet` | *(obrigatório)* | CSV com IDs de amostras e caminhos de arquivos. |
| `--config-file` | *(obrigatório)* | Caminho para o arquivo de configuração YAML a ser gerado. |
| `--output` | *(obrigatório)* | Diretório base de saída. |
| `--run-name` | `undefined` | Nome do run; saída em `{output}/{run_name}/`. |
| `--kraken2-database` | `NA` | Caminho para o banco de dados Kraken2. |
| `--krona-database` | `NA` | Caminho para o banco de taxonomia do Krona. |
| `--taxdump` | `NA` | Diretório do taxdump NCBI (`nodes.dmp`, `names.dmp`). |
| `--host-reference` | `NA` | FASTA do genoma do hospedeiro para remoção com minimap2. |
| `--deacon-index` | `NA` | Índice Deacon para remoção de reads do hospedeiro. |
| `--run-denovo-assembly` | off | Executar montagem de novo com MEGAHIT. |
| `--run-kraken2-reads/--no-kraken2-reads` | on | Kraken2 nas reads. |
| `--run-kraken2-contigs/--no-kraken2-contigs` | on | Kraken2 nos contigs. |
| `--run-diamond-reads`/`--no-diamond-reads` | off | Diamond nas reads. |
| `--run-diamond-contigs`/`--no-diamond-contigs` | off | Diamond nos contigs. |
| `--diamond-database` | `NA` | FASTA de proteínas para Diamond. |
| `--taxids` | `NA` | Mapeamento de taxid NCBI para taxonomia Diamond. |
| `--diamond-sensitivity` | `sensitive` | Modo de sensibilidade do Diamond (`sensitive` / `mid-sensitive` / `more-sensitive` / `ultra-sensitive`). |
| `--evalue` | `0.001` | Limiar de E-value do Diamond. |
| `--bleed-fraction` | `0.005` | Fração do filtro de bleed (max-RPM). |
| `--negative-controls` | (vazio) | IDs de amostras separados por vírgula usados como controles negativos. |
| `--negative-p-threshold` | `0.01` | Limiar de p-valor para filtro de controle negativo. |
| `--minimum-hit-group` | `4` | Parâmetro minimum-hit-group do Kraken2. |
| `--run-reference-assembly`/`--no-run-reference-assembly` | off | Habilitar montagem de genoma baseada em referência ("dynamic reference assembly"). |
| `--method` | `kraken2` | Método alvo para classificar a família viral (`kraken2`, `diamond`, `both`). Obrigatório quando `--run-reference-assembly` está ativo. |
| `--source` | `reads` | Tipo de origem contendo a taxonomia da família (`reads`, `contigs`, `both`). Obrigatório quando `--run-reference-assembly` está ativo. |
| `--reads-count` | `100` | Mínimo de reads mapeadas da família-alvo. |
| `--contigs-count` | `1` | Mínimo de contigs pertencentes a família-alvo. |
| `--families` | `Coronaviridae,...` | Famílias alvo por default, ex: `Coronaviridae,Orthomyxoviridae` |
| `--reference-selection-strategy` | `taxid` | Abordagem para selecionar genoma referência (`taxid` ou `similarity`). |
| `--blast-qcov` | `80` | Cobertura mínima (`--similarity`). |
| `--blast-pident` | `80` | Similaridade mínima (`--similarity`). |
| `--viral-genomes` | `NA` | Caminho do arquivo viral.genomes.fasta. |
| `--viral-taxids` | `NA` | Caminho para mapeamento accession para ID NCBI. |
| `--threads` | `1` | Threads por tarefa. |
| `--threads-total` | `1` | Total de threads para o workflow. |
| `--create-config-only` | off | Apenas gerar o config; não executar o workflow. |

### Exemplos

**Illumina — Kraken2 nas reads (padrão):**

```bash
viralunity meta illumina \
    --sample-sheet /caminho/amostras.csv \
    --config-file /caminho/config_meta.yml \
    --run-name exemplo \
    --output /caminho/saida \
    --kraken2-database /caminho/kraken2_db \
    --krona-database /caminho/krona_taxonomy \
    --taxdump /caminho/taxdump \
    --threads 2 \
    --threads-total 4
```

**Illumina — com adaptadores e remoção do hospedeiro:**

```bash
viralunity meta illumina \
    --sample-sheet amostras.csv \
    --config-file config.yaml \
    --output /caminho/saida \
    --kraken2-database /caminho/kraken2_db \
    --krona-database /caminho/krona_taxonomy \
    --taxdump /caminho/taxdump \
    --adapters /caminho/adaptadores.fa \
    --host-reference /caminho/genoma_hospedeiro.fa \
    --threads 4 --threads-total 8
```

**Illumina — pipeline completo (montagem + Kraken2 + Diamond, reads e contigs):**

```bash
viralunity meta illumina \
    --sample-sheet amostras.csv \
    --config-file config.yaml \
    --output /caminho/saida \
    --kraken2-database /caminho/kraken2_db \
    --krona-database /caminho/krona_taxonomy \
    --taxdump /caminho/taxdump \
    --run-diamond-reads --run-diamond-contigs \
    --diamond-database /caminho/proteinas.faa \
    --taxids /caminho/protein2taxid.tsv \
    --run-denovo-assembly \
    --host-reference /caminho/hospedeiro.fa \
    --threads 4 --threads-total 8
```

**Illumina — com montagem de referência dinâmica:**

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

**Nanopore — pipeline completo:**

```bash
viralunity meta nanopore \
    --sample-sheet /caminho/amostras_nano.csv \
    --config-file /caminho/config_nano.yml \
    --run-name exemplo \
    --output /caminho/saida_nano \
    --kraken2-database /caminho/kraken2_db \
    --krona-database /caminho/krona_taxonomy \
    --taxdump /caminho/taxdump \
    --run-denovo-assembly \
    --run-polish-medaka --medaka-model r941_min_high_g360 \
    --threads 4 --threads-total 4
```

---

## `viralunity get-databases`

O comando `get-databases` baixa e configura os bancos de dados de referência necessários pelo pipeline meta. Cada subcomando recebe um argumento opcional `--path` (padrão: `databases` no diretório atual), dentro do qual um subdiretório nomeado é criado.

```bash
viralunity get-databases --help
viralunity get-databases kraken2 --help
```

### `kraken2` — Índice pré-construído do Kraken2

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--path` | `databases` | Diretório pai; cria `{path}/kraken2/`. |
| `--url` | k2_viral_20240112 | URL do arquivo do índice pré-construído. Verifique os [índices disponíveis](https://benlangmead.github.io/aws-indexes/k2). |

```bash
viralunity get-databases kraken2 --path /dados/dbs
# banco de dados em /dados/dbs/kraken2/
# uso: --kraken2-database /dados/dbs/kraken2
```

### `krona` — Taxonomia do Krona

Requer que o ambiente conda viralunity esteja ativo (`CONDA_PREFIX` deve estar definido).

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--path` | `databases` | Diretório pai; cria `{path}/krona/taxonomy/`. |

```bash
conda activate viralunity
viralunity get-databases krona --path /dados/dbs
# taxonomia em /dados/dbs/krona/taxonomy/
# uso: --krona-database /dados/dbs/krona/taxonomy
```

### `taxdump` — Taxdump NCBI

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--path` | `databases` | Diretório pai; cria `{path}/taxdump/`. |
| `--url` | NCBI taxdump | URL do arquivo taxdump. |

```bash
viralunity get-databases taxdump --path /dados/dbs
# nodes.dmp, names.dmp, etc. em /dados/dbs/taxdump/
# uso: --taxdump /dados/dbs/taxdump
```

### `diamond` — Banco de dados de proteínas Diamond

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--path` | `databases` | Diretório pai; cria `{path}/diamond/`. |
| `--taxon` | Viruses | Nome do táxon NCBI para baixar (ex. `Viruses`, `coronaviridae`). |
| `--refseq`/`--no-refseq` | `on` | Limitar apenas a genomas RefSeq. |
| `--threads` | `1` | Threads para `diamond makedb`. |
| `--skip-makedb` | off | Apenas baixar e reformatar arquivos; não executar `diamond makedb`. |

Requer a CLI NCBI Datasets (`conda install -c conda-forge ncbi-datasets-cli`).

```bash
viralunity get-databases diamond --path /dados/dbs --threads 4
# /dados/dbs/diamond/viral.dmnd
# /dados/dbs/diamond/protein2taxid.tsv
# uso: --diamond-database /dados/dbs/diamond/viral.dmnd
#      --taxids /dados/dbs/diamond/protein2taxid.tsv
```

### `virus-genome` — Banco de dados de genomas virais

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--path` | `databases` | Diretório pai; cria `{path}/virus_genomes/`. |
| `--taxon` | Viruses | Nome do táxon NCBI para baixar (ex. `Viruses`, `coronaviridae`). |
| `--refseq`/`--no-refseq` | `on` | Limitar apenas a genomas RefSeq. |
| `--skip-makeblastdb` | off | Apenas baixar e reformatar arquivos; não executar `makeblastdb`. |

Requer a CLI NCBI Datasets e BLAST+ (`makeblastdb`).

Após baixar e reformatar as sequências, o comando executa automaticamente `makeblastdb` para criar um índice BLAST nucleotídico ao lado do arquivo FASTA. Esse índice é necessário quando a estratégia de seleção de referência por similaridade (`--reference-selection-strategy similarity`) for utilizada.

```bash
viralunity get-databases virus-genome --taxon Viruses --path /dados/dbs
# FASTA em /dados/dbs/virus_genomes/viral.genomes.fasta
# Índice BLAST em /dados/dbs/virus_genomes/viral.genomes.fasta.{nhr,nin,nsq,...}
# taxids em /dados/dbs/virus_genomes/genome2taxid.tsv
# uso: --viral-genomes /dados/dbs/virus_genomes/viral.genomes.fasta
#      --viral-taxids /dados/dbs/virus_genomes/genome2taxid.tsv
```

### `host-genome` — Download de genoma do hospedeiro

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--path` | `databases` | Diretório pai; cria `{path}/host_genomes/`. |
| `--accession` | *(obrigatório)* | ID de acesso NCBI do genoma (ex. `GCA_000001405.29`). |

```bash
viralunity get-databases host-genome --accession GCA_000001405.29
# fasta em databases/host_genomes/GCA_000001405.29.fasta
```

### `deacon-index` — Índice pre-construído do Deacon

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--path` | `databases` | Diretório pai; cria `{path}/deacon_indexes/`. |
| `--index-name` | `panhuman-1` | Nome do índice pré-construído (`panhuman-1` ou `panmouse-1`). |

```bash
viralunity get-databases deacon-index --index-name panhuman-1
# index em databases/deacon_indexes/panhuman-1.idx
# uso: --deacon-index databases/deacon_indexes/panhuman-1.idx
```

### Baixar todos os bancos de dados de uma vez

```bash
viralunity get-databases all --threads 4
```

---

## `viralunity build-deacon-index`

Criar um índice Deacon a partir de um FASTA de genoma do hospedeiro.

| Opção | Padrão | Descrição |
|-------|--------|-----------|
| `--path` | `databases` | Diretório pai; cria `{path}/deacon_indexes/`. |
| `--input` | *(obrigatório)* | Arquivo FASTA de entrada para criar o índice. |
| `--threads` | `8` | Número de threads a utilizar. |

```bash
viralunity build-deacon-index \
    --input databases/host_genomes/GCA_000001405.29.fasta \
    --threads 4
# saída em databases/deacon_indexes/GCA_000001405.29.idx
```

---

## Sobrescrita pelo arquivo de configuração

Alguns parâmetros de baixo nível são configuráveis somente pelo arquivo YAML gerado via `--config-file` (não estão expostos como flags da CLI porque raramente precisam ser alterados). Os valores padrão preservam o comportamento histórico — a maioria dos usuários pode ignorar esta seção.

Abra o YAML gerado após executar com `--create-config-only` e edite a chave correspondente na seção `# parameters`:

| Chave | Padrão | Efeito |
|-------|--------|--------|
| `minimap2_consensus_align_flags` | `-a --sam-hit-only --secondary=no --score-N=0` | Flags passadas ao `minimap2` ao realinhar o consenso por amostra contra a referência para gerar o alinhamento múltiplo final. Usada pelos pipelines `consensus`. |
| `diamond_max_target_seqs` | `1` | Valor de `--max-target-seqs` do DIAMOND. Aumente se precisar de múltiplos hits proteicos por query (deixa a execução mais lenta e aumenta o tamanho da saída). Usada pelo `meta`. |
| `kraken2_extra_flags` | `--report-minimizer-data` | Flags adicionais incluídas em toda chamada do Kraken2, junto com `--threads` e `--minimum-hit-group`. Use `""` para remover a coluna de dados de minimizers do relatório. Usada pelo `meta`. |

Exemplo: para aumentar o DIAMOND para até cinco hits por query e remover a coluna de minimizers do Kraken2, edite o YAML para:

```yaml
diamond_max_target_seqs: 5
kraken2_extra_flags: ""
```

Em seguida, reexecute o mesmo `viralunity meta ...` (sem `--create-config-only`) passando o config editado para o Snakemake.
