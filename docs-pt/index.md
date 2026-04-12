# Documentação ViralUnity

```{admonition} 🌐 Language / Idioma
:class: tip

This documentation is also available in **English**: [Click here](https://viralunity.readthedocs.io/en/latest/)
```

**ViralUnity** é uma ferramenta para análise de dados de sequenciamento de alto desempenho de vírus. Reúne scripts Python e workflows Snakemake, com diversas dependências, para realizar controle de qualidade, classificação taxonômica e montagem de genomas de referência.

O ViralUnity roda em sistemas *nix e é capaz de processar runs de sequenciamento completos em tempo mínimo em um computador comum.

## Funcionalidades principais

- **Montagem de consenso baseada em referência** para dados Illumina e Nanopore
- **Pipeline de metagenômica** com Kraken2, Diamond e montagem de novo (MEGAHIT)
- **Suporte a vírus segmentados** — monte genomas multi-segmento em uma única execução
- **Remoção de reads do hospedeiro** via minimap2 ou Deacon
- **Gerenciamento de bancos de dados** com comandos de download integrados

## Conteúdo da documentação

```{toctree}
:maxdepth: 2

instalacao
uso
comandos
saida
notas
citacao
```

## Links rápidos

- [Repositório GitHub](https://github.com/InstitutoTodosPelaSaude/ViralUnity)
- [Issues / Bugs](https://github.com/InstitutoTodosPelaSaude/ViralUnity/issues)
- **Licença:** MIT
