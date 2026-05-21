# Visão Geral de Uso

O ViralUnity oferece cinco subcomandos principais:

```
viralunity [--version]
├── build-deacon-index   Criar índice Deacon a partir de um arquivo FASTA
├── create-samplesheet   Gerar planilha de amostras a partir de um diretório de run
├── get-databases        Baixar e configurar bancos de dados de referência
│   ├── kraken2          Baixar índice pré-construído do Kraken2
│   ├── krona            Configurar banco de dados de taxonomia do Krona
│   ├── taxdump          Baixar o taxdump do NCBI
│   ├── diamond          Baixar proteínas virais do RefSeq e construir DB Diamond
│   ├── host-genome      Baixar genoma do hospedeiro via NCBI Datasets
│   └── deacon-index     Baixar índice Deacon pré-construído
├── consensus
│   ├── illumina          Montagem de consenso baseada em referência para dados Illumina
│   └── nanopore          Montagem de consenso baseada em referência para dados Nanopore
└── meta
    ├── illumina          Pipeline de metagenômica para dados Illumina
    └── nanopore          Pipeline de metagenômica para dados Nanopore
```

Use `--help` em qualquer nível para ver a lista completa de opções:

```bash
viralunity --help
viralunity consensus --help
viralunity consensus illumina --help
viralunity meta nanopore --help
```

## Fluxo geral

1. **Gere a planilha de amostras** com `viralunity create-samplesheet`
2. **Baixe bancos de dados** (pipeline meta) com `viralunity get-databases`
3. **Execute o pipeline** com `viralunity consensus` ou `viralunity meta`

```{tip}
Sempre utilize caminhos absolutos para evitar erros ao especificar arquivos.
```
