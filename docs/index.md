# ViralUnity Documentation

```{admonition} 🌐 Language / Idioma
:class: tip

Esta documentação também está disponível em **Português**: [Clique aqui](https://viralunity.readthedocs.io/pt-br/latest/)
```

**ViralUnity** is a simple tool to perform analysis of viral high-throughput sequencing data. It comprises a collection of Python scripts and Snakemake workflows, including several software dependencies, to perform data quality control, taxonomic assignments, and reference genome assembly.

ViralUnity runs on *nix systems and is able to process entire sequencing runs in minimal time on a regular computer.

## Main Features

- **Reference-based consensus assembly** for Illumina and Nanopore data
- **Metagenomics** pipeline with Kraken2, Diamond, and de novo assembly (MEGAHIT)
- **Segmented virus support** — assemble multi-segment genomes in a single run
- **Host depletion** via minimap2 or Deacon
- **Flexible database management** with built-in download commands

## Documentation Contents

```{toctree}
:maxdepth: 2

installation
usage
commands
output
notes
citation
```

## Quick Links

- [GitHub Repository](https://github.com/InstitutoTodosPelaSaude/ViralUnity)
- [Issues / Bugs](https://github.com/InstitutoTodosPelaSaude/ViralUnity/issues)
- **License:** MIT
