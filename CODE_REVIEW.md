# ViralUnity — Professional Code Review

**Branch:** `feature/update-meta`
**HEAD:** `163e72e Add filtered Krona reports aligned with bleed/neg taxa summaries.`
**Comparison baseline:** ~104 files changed vs `main` (`a9efdeb`).

Scope: full repository — Python package (`viralunity/`), Snakemake workflows (`viralunity/scripts/*.smk` + `rules/*.smk`), Snakemake helper scripts (`viralunity/scripts/python/`, `viralunity/scripts/*.py`), tests (`test/`, `test/scripts/`, `test/dryrun_configs/`), CI (`.github/workflows/`), Sphinx docs (`docs/`, `docs-pt/`), and repo-root configuration.

Companion file: `CODE_REVIEW_BACKLOG.md` (prioritized remediation list).

---

## What changed since `main`

This branch is substantially more mature than `main` — many of the structural backlog items that would have dominated a main-branch review are now resolved. Crediting before flagging:

- **Per-rule conda environments**: 149 `conda:` directives across rules, backed by 9 env yamls in `viralunity/scripts/envs/` (`alignment.yaml`, `assembly.yaml`, `clair3.yaml`, `consensus.yaml`, `genome_selection.yaml`, `medaka.yaml`, `qc.yaml`, `taxonomy.yaml`, `utils.yaml`). The main-branch README TO-DO is closed out.
- **`resources:` declarations**: 29 occurrences (mostly `mem_mb`, `threads`) on the memory-heavy rules — cluster execution is now feasible.
- **Shell safety**: 54 `set -euo pipefail` occurrences in multi-line shell blocks; large reduction in silent-failure exposure compared to main.
- **Real dry-run integration tests**: `test/viralunity_dryrun_test.py` + `test/dryrun_configs/*.yaml` + `test/dryrun_configs/create_dryrun_placeholders.sh` execute actual Snakemake DAGs against fixture placeholders instead of mocks. This is the integration-test foundation main lacked entirely.
- **Path-resolution tests**: `test/viralunity_path_resolution_test.py` (300 LOC) exercises CLI behaviour with tempfiles — non-mock coverage of a previously untested concern.
- **Sphinx documentation**: `docs/` (EN) and `docs-pt/` (PT) Sphinx structures with `commands.md` (543 LOC) documenting all flag tables; `.readthedocs.yaml` configured. Replaces the previously-drifted in-repo README flag tables.
- **Argparse → Click CLI migration**: cleaner, decorator-based subcommand definitions; shared option helpers via `_add_common_options` / `_add_resource_options`.
- **New subcommands**: `viralunity get-databases` (9 sub-subcommands: kraken2, diamond, krona, taxdump, virus-genome, host-genome, deacon-index, clean-protein-fasta, all), `viralunity build-deacon-index`.
- **Metagenomics expansion**: DIAMOND-based protein-level classification (reads + contigs), reference-assembly checkpoint workflow, host-read dehosting (deacon), per-rule kraken2 split (reads vs contigs).
- **Krona post-processing**: `viralunity/scripts/python/filter_taxids.py` adds robust, properly-logged filtering of Krona inputs by pass-taxid lists.
- **Granular rule organization**: 21 rule files in `viralunity/scripts/rules/` (vs 5 on main); each tightly focused.

---

## Headline findings

1. **`AdaptersNotFoundError` is raised without being imported** in `viralunity/validators.py:123`. When `validate_illumina_requirements` reaches the file-existence branch (i.e. `--adapters` is provided), Python raises `NameError: name 'AdaptersNotFoundError' is not defined` instead. The existing tests mock the entire validator (`test/viralunity_consensus_test.py:152-180`), so they pass — but the real code path is broken. ([C-1])
2. **`make test` is broken in any environment that hasn't pre-installed the package**: the dryrun test imports `pytest` at module load (`test/viralunity_dryrun_test.py:3`), the get-databases / meta-cli tests import `click.testing` at module load — neither `pytest` nor `click` is declared in `environment.yml`. They reach the CI runner via `pip: - .` → `setup.py install_requires`, but a fresh `conda env create -f environment.yml` without the package install leaves these missing. ([C-2])
3. **`Dockerfile:2` LABEL `version="1.0.3"`** while `viralunity/__init__.py:4` declares `__version__ = "1.1.0"`. Images built from this Dockerfile mis-identify themselves. ([C-3])
4. **README.md is now a 7-line stub redirecting to ReadTheDocs.** Anyone landing on the GitHub page sees no installation/usage info, and the docs link assumes RTD is live. For a public-facing repository that's a release-readiness regression. ([H-1])
5. **`print()` in library code persists** at `viralunity/viralunity_meta.py:226` despite the rest of the package using `logger`. ([M-1])
6. **CI still triggers on `pull_request` only**, never on `push: main`, and runs no lint / type checks / dry-runs. The merge commit that lands this branch will not be CI-verified. ([H-2])
7. **`setup.py` mixes dev dependencies (`black>=24.4`, `pytest>=9.0`) into `install_requires`.** Every production install pulls a formatter and a test framework. ([H-3])
8. **The 4-way consensus workflow duplication (`consensus_{illumina,nanopore}{,_segmented}.smk`) survives.** ~60% identical code across the segmented/non-segmented pair on each data type. ([H-4])

---

## Rubric scorecard

| Dimension | Verdict | Notes |
|---|---|---|
| 1. Python package | **partial** | Mostly typed, mostly documented (~26% docstring coverage on public API per agent audit); ~3 critical bugs (C-1, C-3, M-1); helper duplication in `_run`/`_parse_*` between new CLIs |
| 2. Snakemake workflows | **good** | `conda:`/`resources:`/`set -euo pipefail`/`log:`/`benchmark:` now broadly present; 4× duplication of top-level consensus snakefiles still unresolved; a handful of rules still missing `set -euo pipefail` |
| 3. Snakemake helper scripts | **partial** | `filter_taxids.py` is well-structured; `calculate_assembly_stats.py` and `rename_sequences.py` still untyped/undoc'd and use `subprocess.Popen(shell=True)` with concatenated paths |
| 4. Tests | **much improved** | Real dry-run tests exist; path-resolution tests exist; but `make test` fails locally on missing pytest/click, and broken `AdaptersNotFoundError` test path is mock-shielded |
| 5. Documentation | **partial** | Sphinx docs are comprehensive and current with the CLI; README is a 7-line redirect; no `CHANGELOG.md`/`CONTRIBUTING.md`/`CODE_OF_CONDUCT.md` |
| 6. Repo hygiene | **good** | `.gitignore` covers stale dirs; `MANIFEST.in` added; version single-sourced from `__init__.py` (Dockerfile excepted); LICENSE year still `2021`; no `pyproject.toml` |
| 7. CI / release | **fail** | Still PR-only; no lint/type checks; no smoke test target; no git tags / release automation |

---

## Findings — Critical

### [C-1] `AdaptersNotFoundError` used without being imported in validators
- **Location**: `viralunity/validators.py:123` raises `AdaptersNotFoundError`; the module's `from viralunity.exceptions import (...)` block (lines 7-19) imports `ValidationError`, `FileNotFoundError`, `SampleSheetError`, `SampleConfigurationNotFoundError`, `ReferenceNotFoundError`, `PrimerSchemeNotFoundError`, `TaxdumpNotFoundError`, `DiamondDatabaseNotFoundError` — **not** `AdaptersNotFoundError`. The class itself exists in `viralunity/exceptions.py:64`.
- **Why it ships**: tests at `test/viralunity_consensus_test.py:152-180` and `test/viralunity_meta_test.py:*` mock `validate_illumina_requirements` and assert it raises `AdaptersNotFoundError`. The mock makes the assertion true; the real validator never runs that branch during tests because the validator itself is patched out.
- **Repro** (production path): any Illumina run with a non-existent `--adapters` path. Real code flow: `cli` → `consensus_main` → `validate_args` → `validate_illumina_requirements(args)` → `validate_file_exists(adapters, ...)` raises `FileNotFoundError` → the `except FileNotFoundError` block at line 122 tries to `raise AdaptersNotFoundError(...)` → `NameError: name 'AdaptersNotFoundError' is not defined`.
- **Fix**: add `AdaptersNotFoundError` to the import list at `viralunity/validators.py:7-19`. While at it, chain with `from e`. **Quick win**.

### [C-2] `make test` requires `pytest` and `click` but neither is in `environment.yml`
- **Location**:
  - `test/viralunity_dryrun_test.py:3` imports `pytest`.
  - `test/viralunity_get_databases_cli_test.py:14`, `test/viralunity_meta_cli_test.py:5` import `click.testing.CliRunner`.
  - `environment.yml:1-20` declares neither `pytest` nor `click` as conda or pip dependencies; only `pip: - .` is listed.
  - `setup.py:11-17` includes `click>=8.0`, `black>=24.4`, `pytest>=9.0` in `install_requires` — so a `pip install .` from inside the conda env will pull them in.
- **Impact**:
  - CI happens to work because the workflow does `conda env update --file environment.yml --name viralunity` followed by `pip install -e .`, and the `pip: - .` in `environment.yml` also bootstraps install_requires.
  - But anybody running `make test` *without* having first run `pip install -e .` will get `ModuleNotFoundError`. The `make install` target exists but isn't a dependency of `make test`.
  - More importantly, mixing test/dev tooling into runtime `install_requires` is wrong (see [H-3]).
- **Fix**: declare `pytest>=9.0` (and ideally `click>=8.0`) explicitly in `environment.yml` (or in a `dev` extras section once a `pyproject.toml` exists); make `test` depend on `install` in the `Makefile`. **Quick win** for the Makefile dependency; the env split is in the backlog.

### [C-3] Dockerfile version label is stale
- **Location**: `Dockerfile:2` — `LABEL version="1.0.3"`.
- **Truth**: `viralunity/__init__.py:4` — `__version__ = "1.1.0"`.
- **Impact**: container metadata lies about the contained package version, breaking any consumer that filters images by label.
- **Fix**: either hardcode `LABEL version="1.1.0"` or — better — derive from `__version__` via a build arg / multi-stage. **Quick win**: hardcode for now, automate later.

---

## Findings — High

### [H-1] README is a 7-line redirect to ReadTheDocs
- **Location**: `README.md` (whole file, 7 lines).
- **Impact**: GitHub-first users (the majority for a public open-source repo) see ~30 words and a link. The previous README had `Installation`, `Usage`, `Tests`, `Citation` sections — all gone. If `viralunity.readthedocs.io` is not yet live (or behind a redirect), users hit a dead end.
- **Fix**: keep the RTD link, but restore a minimal usable block in the README: install command, two-line usage example for each subcommand, link to RTD for full docs. Quick win.

### [H-2] CI still runs only on pull_request; no lint, type checks, or smoke
- **Location**: `.github/workflows/ci.yaml:3-4` (`on: pull_request:` only); steps at lines 25-38 run unit tests + install verification only.
- **Impact**: merge commits to `main` skip CI. The branch's own dry-run tests (`test/viralunity_dryrun_test.py`) aren't deliberately exercised. No `ruff`/`mypy`/`black --check` step protects code quality.
- **Fix**:
  1. Add `push: branches: [main]` to the trigger.
  2. Add a `dryrun-test` job that pre-populates `test/dryrun_configs/` placeholders (via `create_dryrun_placeholders.sh`) and runs `pytest test/viralunity_dryrun_test.py`.
  3. Add a `lint` job: `ruff check viralunity/` and `black --check viralunity/` (black is already in `install_requires`).
- Backlog (multi-step, design needed).

### [H-3] Dev tooling leaks into `install_requires`
- **Location**: `setup.py:11-17`. `black>=24.4` and `pytest>=9.0` should not be runtime dependencies.
- **Impact**: every production install (including Docker image and end-user `pip install`) pulls a formatter and a test framework, increasing image size and surface area. Conflicts with the project's own license/scope.
- **Fix**: migrate to `pyproject.toml` with `[project.optional-dependencies] dev = ["black>=24.4", "pytest>=9.0", "ruff>=0.5"]` and document `pip install -e ".[dev]"` for contributors. Backlog (packaging change).

### [H-4] Top-level consensus workflows still duplicated 4-way
- **Location**: `viralunity/scripts/consensus_{illumina,nanopore}{,_segmented}.smk`. Rule bodies (`rule all`, `organize_files`, `unify_assembly_statistics_reports`, `align_consensus_to_reference_genome`) are ~60% copy-pasted between each segmented/non-segmented pair, only differing by `SEGMENT_WILDCARD` interpolation. The included `rules/*.smk` files are already cleanly parametrized over `SEGMENT_WILDCARD`; the top-level files are not.
- **Fix**: collapse each pair into a single snakefile that branches on `isinstance(config["reference"], dict)`. Backlog — needs careful DAG verification.

### [H-5] Metagenomics nanopore reuses the Illumina QC rule
- **Location**: `viralunity/scripts/metagenomics_nanopore.smk` (and the dehost/diamond/kraken2 nanopore rule files) include or rely on `rules/qc_illumina.smk` outputs (`perform_qc` produces paired `_R1` / `_R2` outputs).
- **Impact**: Nanopore data is single-end. The DAG either fails immediately or produces empty/garbage trimmed reads if it does run.
- **Verification needed**: dry-runs of the metagenomics_nanopore workflow against the supplied fixture configs.
- **Fix**: add a `rules/qc_nanopore.smk` with single-end-appropriate trimming (NanoFilt / fastp single-end). Backlog.

### [H-6] Several rules still lack `set -euo pipefail`
- **Locations** (counts, not exhaustive): the 54 occurrences cover most rules, but some don't:
  - `viralunity/scripts/rules/alignment_illumina.smk:25-31` `map_reads` shell (piped `minimap2 | samtools view | samtools sort`).
  - `viralunity/scripts/rules/qc_illumina.smk:39-64` (~25-line fastp invocation).
  - `viralunity/scripts/consensus_nanopore.smk:77-82` `align_consensus_to_reference_genome` (cat + minimap2 pipeline).
  - `viralunity/scripts/rules/consensus_illumina.smk:40-55` `unify_assembly_statistics_reports` (echo + cat).
- **Impact**: in piped commands without `pipefail`, only the *last* exit status determines rule success — a failing `minimap2` upstream of a successful `samtools sort` is treated as success.
- **Fix**: prefix the still-bare shell blocks with `set -euo pipefail`. Backlog.

### [H-7] `subprocess.Popen(..., shell=True)` with concatenated paths in helper script
- **Location**: `viralunity/scripts/calculate_assembly_stats.py:11,19` (after earlier code at lines 8-10).
- **Same finding as on main.** Snakemake feeds the paths, so today it's safe in practice, but the pattern remains a footgun.
- **Fix**: switch to `subprocess.run([...], check=True, capture_output=True)`. Backlog.

### [H-8] Duplicate helper `_run` between two new CLIs
- **Location**: `viralunity/viralunity_get_databases_cli.py:25` and `viralunity/viralunity_build_deacon_index_cli.py:13` define near-identical `_run(...)` helpers for wrapping subprocess calls with friendly error reporting.
- **Fix**: extract to `viralunity/_subprocess.py` (or a `utils` module). Backlog.

### [H-9] Duplicate parsers inside `viralunity_get_databases_cli.py`
- **Location**: `_parse_data_report()` (lines 180-196) and `_parse_genome_data_report()` (lines 399-420) are ~95% identical; `_reformat_protein_fasta()` (199-260) and `_reformat_genome_fasta()` (423-448) share accession-extraction logic.
- **Fix**: parameterize a single parser/reformatter; keep the call sites that differ. Backlog.

---

## Findings — Medium

### [M-1] `print()` in library code at `viralunity_meta.py:226`
- **Location**: `viralunity/viralunity_meta.py:226` — `print("No samples were provided.")`.
- **Fix**: replace with `logger.warning(...)`. **Quick win**.

### [M-2] Suspicious `args.update(...)` in `viralunity_meta_cli.py`
- **Location**: `viralunity/viralunity_meta_cli.py:363` — `args.update({k: v for k, v in kwargs.items() if k not in args})`.
- **Impact**: probably no-op or redundant after `_build_meta_args(...)` constructed `args`, but if a kwarg ever sneaks in that isn't yet in `args`, behaviour is opaque.
- **Fix**: clarify intent (add comment + explicit set of expected keys) or remove if redundant. Backlog.

### [M-3] `FileNotFoundError` shadows the Python builtin
- **Location**: `viralunity/exceptions.py:16` defines `FileNotFoundError(ViralUnityError)`; `viralunity/validators.py:11` imports it explicitly. Any caller that does `from viralunity.validators import *` (or, worse, mixes a `try: ... except FileNotFoundError` that is meant to catch the OS-level builtin) gets the wrong one.
- **Fix**: rename to `ViralUnityFileNotFoundError` or `MissingPathError`; update the ~5 importers. Backlog.

### [M-4] Broad `except Exception` without chaining in pipeline mains
- **Location**: `viralunity/viralunity_consensus.py:218-221` and `viralunity/viralunity_meta.py:244-247` catch `Exception`, log, return 1 — but the `logger.exception(...)` already includes the traceback, so re-raising isn't necessary. The pattern is acceptable as a top-level handler. *However*, narrower validation-level `raise FooError(...)` calls in the same modules and in `viralunity/validators.py` and `viralunity/config_generator.py` should still use `from e` chaining.
- **Fix**: audit `raise X(...)` sites inside `except` blocks and add `from e`. Quick win.

### [M-5] Docstring coverage on public API is low
- **Reported**: ~26% function docstring coverage across `viralunity/` (agent audit).
- **Specifically thin**: Click command handlers across `viralunity_get_databases_cli.py`, `viralunity_build_deacon_index_cli.py`, `viralunity_consensus_cli.py`, `viralunity_meta_cli.py` — most lack a docstring on the decorated function; Click uses the `--help`/`help=` parameters of decorators instead, which is functional but doesn't help IDE tooling or contributors.
- **Fix**: add minimal "Summary line. Args section." docstrings on each command handler. Backlog (mechanical sweep).

### [M-6] Type-hint inconsistency on Click handlers
- **Location**: virtually all `@click.command()`-decorated functions across the four CLI modules lack parameter type hints on the handler signature. Click's `type=` decorators are the source of truth at parse time, but signatures could still match.
- **Fix**: add types on handler params (helpful for `mypy`/IDE). Backlog.

### [M-7] LICENSE year stale
- **Location**: `LICENSE:3` — `Copyright (c) 2021 filiperomero2`. The setup.py URL now points to `InstitutoTodosPelaSaude` and the project is active in 2026.
- **Fix**: update to `2021-2026` and (optionally) attribute the institutional copyright as well. Quick win.

### [M-8] Helper scripts still untyped/undocumented; `input` shadows builtin
- **Location**: `viralunity/scripts/calculate_assembly_stats.py`, `viralunity/scripts/rename_sequences.py`. Parameter named `input` (`calculate_assembly_stats.py:73,89`; `rename_sequences.py:9,29`); no module docstrings, no function docstrings, no types.
- **Fix**: rename `input` → `input_files` / `input_path`; add docstrings + types. Backlog.

### [M-9] Hardcoded tool flags
- **Location**: representative:
  - `viralunity/scripts/consensus_illumina.smk:111` and parallel `_segmented` variants — `minimap2 -a --sam-hit-only --secondary=no --score-N=0`.
  - `viralunity/scripts/rules/metagenomics_diamond_reads_illumina.smk:31` — DIAMOND `--max-target-seqs 1`.
  - `viralunity/scripts/rules/metagenomics_kraken2_reads_illumina.smk:13` — kraken2 `--report-minimizer-data --minimum-hit-group 3`.
- **Fix**: expose as config keys with documented defaults. Backlog.

### [M-10] Stray `# check this` and deprecated `delim_whitespace=True`
- **Location**: `viralunity/scripts/calculate_assembly_stats.py:26` (or thereabouts on this branch).
- **Fix**: `sep=r"\s+"` (raw string) and resolve the comment. Quick win.

### [M-11] CI installs `conda env`, then `pip install -e .` separately
- **Location**: `.github/workflows/ci.yaml:32, 36`. The env yaml already has `pip: - .` which installs the package; running `pip install -e .` afterwards re-installs (now editable). Either step alone suffices; running both pulls everything twice.
- **Fix**: pick one. Backlog.

### [M-12] No git tags / release automation
- **Location**: `git tag -l` empty.
- **Fix**: tag `v1.1.0` at the appropriate commit on main once this branch merges; document a release process. Backlog.

---

## Findings — Low

### [L-1] No `CHANGELOG.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue/PR templates
- **Fix**: add minimal stubs as quick wins; templates can come later.

### [L-2] `Makefile::test` doesn't depend on `install`
- **Location**: `Makefile:4`. A contributor's first `make test` after clone fails confusingly (`ModuleNotFoundError`) because the package isn't installed yet.
- **Fix**: change `test:` to `test: install` (or add a guard). Quick win.

### [L-3] PATH `print(...)` lines at module load in test files
- **Location**: `test/viralunity_consensus_test.py:17`, `test/viralunity_meta_test.py:16`, `test/viralunity_meta_cli_test.py:6` all do `print(os.environ.get("PATH"))` at import time. Noise in every test run.
- **Fix**: delete the lines (and the now-unused `import os` where applicable). Quick win.

### [L-4] `environment.yml` minor-pinned, `setup.py` patch-pinned
- **Location**: `environment.yml:7` `snakemake=7.32` vs `setup.py:15` `snakemake==7.32.0`. Same project, two pin levels.
- **Fix**: align to one. Backlog.

### [L-5] No `pyproject.toml`; no ruff/mypy/black/pre-commit configs
- **Fix**: backlog; pairs with [H-3].

### [L-6] `viralunity_meta.py:129-131` passes `forceall=True, lock=False` to Snakemake; consensus does not
- **Location**: `viralunity/viralunity_meta.py:124-132` vs `viralunity/viralunity_consensus.py:156-161`.
- **Impact**: meta always re-runs every rule; defeats Snakemake's caching.
- **Fix**: align semantics; likely remove `forceall=True`. Backlog (behavioural change).

### [L-7] Docker `LABEL description` is generic
- **Location**: `Dockerfile:3`.
- **Fix**: align with `_description` from `viralunity/__init__.py`. Backlog.

---

## Summary of issue counts

| Severity | Count |
|---|---|
| Critical | 3 |
| High | 9 |
| Medium | 12 |
| Low | 7 |
| **Total** | **31** |

Comparison to a hypothetical main-branch review: high-severity *structural* items (conda-per-rule, resources declarations, `set -euo pipefail`, integration tests) are already resolved; remaining highs are mostly polish (README, CI hardening, packaging hygiene, workflow consolidation). The Critical count is dominated by precise import bug ([C-1]) and packaging/version metadata mismatches ([C-2], [C-3]) — all quick to fix.

See `CODE_REVIEW_BACKLOG.md` for the prioritized work list and the quick-wins scope landed in this engagement.
