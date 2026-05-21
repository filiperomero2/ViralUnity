# Changelog

All notable changes to this project will be documented in this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aspires to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The release process is documented in [RELEASING.md](RELEASING.md).

## [1.1.0] - 2026-05-21

### Added

- New `viralunity get-databases` subcommand group with downloaders for
  kraken2, krona, taxdump, virus-genome, host-genome, deacon-index, and
  diamond databases, plus an `all` aggregator and a `clean-protein-fasta`
  utility.
- New `viralunity build-deacon-index` subcommand to build a Deacon
  minimizer index from a host FASTA.
- DIAMOND-based protein-level classification (reads + contigs) for both
  Illumina and Nanopore metagenomics workflows.
- Reference-assembly checkpoint workflow that builds per-family
  consensuses from filtered Kraken2/DIAMOND hits.
- Filtered Krona reports aligned with bleed and negative-control taxa
  summaries.
- Per-rule conda environments (`viralunity/scripts/envs/*.yaml`) and
  `resources:` declarations on memory-heavy rules; cluster execution is
  now feasible.
- Real Snakemake dry-run integration tests
  (`test/viralunity_dryrun_test.py` + `test/dryrun_configs/`).
- Path-resolution unit tests (`test/viralunity_path_resolution_test.py`).
- Sphinx documentation under `docs/` (EN) and `docs-pt/` (PT) with full
  CLI flag tables; ReadTheDocs configuration.
- `pyproject.toml` migration (PEP 621): `[project]` metadata, dev extras
  (`black`, `pytest`, `ruff`), `[tool.black]`, `[tool.ruff]`,
  `[tool.mypy]`, `[tool.pytest.ini_options]` configurations.
- CI hardening: `push: main` trigger, dedicated `lint` job (black --check,
  ruff check), dedicated `dryrun` job; CI now runs on PRs *and* on
  pushes to `main`.
- Three previously-hardcoded tool flags are now config-driven keys
  (`minimap2_consensus_align_flags`, `diamond_max_target_seqs`,
  `kraken2_extra_flags`) with the historical defaults preserved.
- Community files: `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` (Contributor
  Covenant 2.1), `.github/ISSUE_TEMPLATE/`, `.github/PULL_REQUEST_TEMPLATE.md`.
- `RELEASING.md` documenting the version-bump and tag workflow.

### Changed

- Migrated from `argparse` to Click for all CLI surfaces; shared options
  are stacked via `_add_common_options` and `_add_resource_options`
  decorators.
- Granular Snakemake rule organisation: 21 files under
  `viralunity/scripts/rules/` (vs 5 on `main`).
- Aligned dependency pinning between `environment.yml` and
  `pyproject.toml` on `>=X.Y` semantics, with `<8` upper bound on
  Snakemake (its 8.x release dropped the Python API this package calls).
- LICENSE copyright extended to `2021-2026`.
- README restored to a usable shape: install + quick-start blocks for
  each subcommand, with the ReadTheDocs link kept as the canonical
  source.
- Dockerfile `LABEL description` now matches `viralunity._description`
  ("A pipeline for viral metagenomics analysis.") instead of the
  generic "Docker image for viralunity" stub.
- `viralunity meta` no longer passes `forceall=True, lock=False,
  workdir=os.getcwd()` to `snakemake(...)`. The pipeline now respects
  Snakemake's caching and matches the consensus pipeline's invocation;
  the workdir defaults to cwd, which is what `resolve_path_args`
  already prepares paths against.

### Fixed

- **Critical**: `viralunity/validators.py` raised `AdaptersNotFoundError`
  without importing it; any Illumina run providing a non-existent
  `--adapters` path triggered a `NameError` instead of the intended
  `AdaptersNotFoundError`. Existing tests mocked the validator so this
  did not surface in CI.
- Dockerfile `LABEL version` bumped from `1.0.3` to `1.1.0` to match
  `viralunity/__init__.py:__version__`.
- Replaced `print("No samples were provided.")` in
  `viralunity/viralunity_meta.py` with `logger.warning(...)`.
- Added `from e` exception chaining to all custom-exception re-raises
  inside `except ... as e` blocks in `viralunity/validators.py` and
  `viralunity/config_generator.py`.
- Replaced deprecated `pd.read_csv(..., delim_whitespace=True)` with
  `sep=r"\s+"` in
  `viralunity/scripts/python/calculate_assembly_stats.py`.
- Removed `subprocess.Popen(..., shell=True)` from
  `calculate_assembly_stats.py`; now uses `subprocess.run([...])` with
  list args so the path is not interpreted by a shell.
- Added `set -euo pipefail` to every multi-line shell block across the
  rule files so piped commands fail fast instead of letting an upstream
  error get swallowed by a downstream `samtools sort`.
- `make test` now depends on `make install` so a fresh clone runs
  cleanly.
- Removed leftover `print(os.environ.get("PATH"))` debug calls from
  `test/viralunity_consensus_test.py` and `test/viralunity_meta_test.py`.

### Refactored

- Renamed the custom `FileNotFoundError` exception to
  `ViralUnityFileNotFoundError` so it no longer shadows Python's
  builtin.
- Extracted a shared `viralunity._subprocess.run_command` helper out of
  the duplicated `_run` helpers in
  `viralunity_get_databases_cli.py` and
  `viralunity_build_deacon_index_cli.py`.
- Consolidated the two near-identical `data_report.jsonl` parsers in
  `viralunity_get_databases_cli.py` into a single parametric function
  driven by a `key_extractor` callable.
- Extracted shared consensus rules (`generate_multiqc_report`,
  `calculate_assembly_statistics`, `align_consensus_to_reference_genome`)
  into `rules/consensus_illumina_common.smk` and
  `rules/consensus_nanopore_common.smk`; the four top-level snakefiles
  now keep only the segment-specific rules.
- Lifted shared pipeline orchestration (`ConfigGenerator` preamble, the
  `snakemake(...)` invocation, the `main` try/except skeleton) into
  `viralunity/_orchestrator.py`.
- Expanded one-line docstrings on all consensus/meta Click handlers and
  added type hints to all 18 Click handler signatures.
- Polished helper scripts (`calculate_assembly_stats.py`,
  `rename_sequences.py`): module + function docstrings, type hints,
  renamed the `input` parameter that shadowed the Python builtin.
- Removed two no-op `args.update({k: v for k, v in kwargs.items() if k
  not in args})` lines in `viralunity_meta_cli.py` that obscured the
  intent of `_build_meta_args`.

## [1.1.0] — branch-vs-main overview

See git history for the substantial changes on this branch versus
`main`: per-rule conda environments, resources declarations, set -euo
pipefail sweep, dryrun integration tests, sphinx docs, argparse→Click
CLI migration, new `get-databases` / `build-deacon-index` subcommands,
metagenomics expansion (DIAMOND, reference assembly, dehosting,
kraken2 reads/contigs).
