# Contributing to ViralUnity

Thanks for your interest in improving ViralUnity. This document covers
the basics of setting up a development environment, running the test
suite, keeping the code lint-clean, and getting changes merged.

## Setup

ViralUnity targets Python 3.8 – 3.11 and is distributed alongside a
conda environment that pins the bioinformatics tooling (snakemake,
minimap2, fastp, kraken2, diamond, krona, deacon, blast, …).

```bash
# Clone and create the conda environment
git clone https://github.com/InstitutoTodosPelaSaude/ViralUnity.git
cd ViralUnity
conda env create -f environment.yml
conda activate viralunity

# Install ViralUnity in editable mode with dev dependencies
pip install -e ".[dev]"

viralunity --version
```

The `dev` extras install `black`, `pytest`, and `ruff` (see
`pyproject.toml`'s `[project.optional-dependencies]`).

## Running tests

```bash
# Full unit suite (unittest discover under the hood)
make test

# Snakemake dry-run integration tests (require the conda env to be active
# so the bundled tools resolve). Six dryrun configs cover
# consensus_{illumina,nanopore}{,_segmented} and
# metagenomics_{illumina,nanopore}.
pytest test/viralunity_dryrun_test.py -v
```

If `pytest` fails inside the dryrun tests with
`AttributeError: module 'pulp' has no attribute 'list_solvers'`, the
local `pulp` is too new for Snakemake 7.32; pin it with
`pip install 'pulp<3'`.

## Linting and formatting

```bash
black viralunity/ test/
ruff check viralunity/ test/
```

CI runs both with `--check` and will fail the PR on diffs. The
canonical configurations live in `pyproject.toml` (`[tool.black]`,
`[tool.ruff]`, `[tool.ruff.lint]`).

## Submitting changes

1. Branch from `main` (or, for in-flight work, the relevant feature
   branch). Keep each commit focused on one logical change so the
   history stays bisectable.
2. Open a pull request against `main`. If your change addresses an
   existing issue, link it in the PR description.
3. Make sure CI is green: lint job (`black --check`, `ruff check`),
   unit tests, and dryrun tests must all pass.
4. A maintainer will review and merge.

## Releasing

Cutting a versioned release (bumping `__version__`, tagging, pushing
the image) is documented in [`RELEASING.md`](RELEASING.md).

## Code of conduct

Participation in this project is governed by our
[Code of Conduct](CODE_OF_CONDUCT.md) (Contributor Covenant 2.1).
