# Installation

ViralUnity is a Python package that launches Snakemake workflows. The tool dependencies are documented in the conda environment file `environment.yml` and the workflow-specific dependencies are in `viralunity/scripts/envs/`.

## Clone and create the environment

```bash
git clone https://github.com/InstitutoTodosPelaSaude/ViralUnity.git
cd ViralUnity/
conda env create -f environment.yml
conda activate viralunity
```

Or with **micromamba** (recommended on macOS with Apple Silicon):

```bash
micromamba env create -f environment.yml --platform osx-64
micromamba activate viralunity
```

```{warning}
On macOS with Apple Silicon (M1 or later), the `viralunity/scripts/envs/clair3.yaml` environment may not install correctly due to compatibility constraints in the clair3 dependencies.
```

## Verify the installation

```bash
viralunity --version
viralunity --help
```

## Development install

To work on ViralUnity itself, install the optional `dev` extras (linters and tests):

```bash
pip install -e ".[dev]"
```

See [CONTRIBUTING.md](https://github.com/InstitutoTodosPelaSaude/ViralUnity/blob/main/CONTRIBUTING.md) for the full development workflow.
