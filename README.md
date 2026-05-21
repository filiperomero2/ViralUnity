# ViralUnity

ViralUnity is a tool for analysing viral high-throughput sequencing data. It is a Python package that orchestrates Snakemake workflows for data quality control, taxonomic assignment, and reference genome assembly. ViralUnity runs on *nix systems and can process entire sequencing runs in minimal time on a regular computer.

> **Full documentation:** <https://viralunity.readthedocs.io/en/latest/>

## Installation

```bash
git clone https://github.com/InstitutoTodosPelaSaude/ViralUnity.git
cd ViralUnity
conda env create -n viralunity -f environment.yml
conda activate viralunity
pip install -e .
```

Per-rule conda environments under `viralunity/scripts/envs/` are managed automatically by Snakemake; the top-level `environment.yml` only installs ViralUnity itself and its core runtime dependencies.

## Quick start

Five top-level subcommands are exposed via the `viralunity` CLI:

```bash
viralunity create-samplesheet --input <runs-dir> --output samples.csv
viralunity get-databases all --output databases/
viralunity consensus illumina --sample-sheet samples.csv --reference ref.fasta --output run/
viralunity meta      illumina --sample-sheet samples.csv --kraken2-database <db> --output run/
viralunity build-deacon-index --input host.fasta --output host.dcn
```

Each subcommand has its own `--help`; the same information is exhaustively documented in the `docs/` Sphinx site (rendered on ReadTheDocs at the link above).

## Tests

```bash
make test
```

This installs the package in editable mode (if not already installed) and runs the `unittest` suite under `test/`. Snakemake dry-run tests live in `test/viralunity_dryrun_test.py` and use `pytest`.

## Citation

A scientific publication describing this pipeline is being prepared. Meanwhile, please cite this repository. Primary references for upstream tools (Trimmomatic, FastQC, MultiQC, Minimap2, Samtools, iVar, BEDtools, Kraken2, Krona, DIAMOND, fastp, Clair3, Medaka, Deacon) are listed in the ReadTheDocs site.

## License

MIT — see `LICENSE`.
