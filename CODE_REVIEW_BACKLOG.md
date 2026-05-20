# ViralUnity — Code Review Backlog

**Branch:** `feature/update-meta` (HEAD `163e72e`)
Prioritized checklist of remediation items distilled from `CODE_REVIEW.md`. Severity first, then effort (small → large) within each severity.

---

## Phase 1 — Quick wins (in scope for this engagement)

Mechanical, low-risk, no design discussion needed. Will land in a single follow-up commit on `feature/update-meta`.

- [ ] **[C-1] Add `AdaptersNotFoundError` to the import list** in `viralunity/validators.py:7-19`. Add `from e` chaining at line 122-125 while we're there.
- [ ] **[C-3] Bump Dockerfile version label** at `Dockerfile:2` from `1.0.3` to `1.1.0`.
- [ ] **[M-1] Replace `print(...)` in library code** at `viralunity/viralunity_meta.py:226` with `logger.warning`.
- [ ] **[M-4] Add `from e` exception chaining** at all `raise X(...)` inside `except ... as e` blocks in `viralunity/validators.py` and `viralunity/config_generator.py`.
- [ ] **[M-7] Bump `LICENSE`** copyright to `2021-2026`.
- [ ] **[M-10] Fix deprecated pandas idiom** in `viralunity/scripts/calculate_assembly_stats.py`: replace `delim_whitespace=True` with `sep=r"\s+"` (and resolve the `# check this` comment by removing it).
- [ ] **[H-1] Restore a minimal usable README**: keep the RTD link, but add 30-50 lines covering install, two-line usage for each subcommand, link to RTD for full docs.
- [ ] **[L-1] Add minimal `CHANGELOG.md`** documenting the quick wins.
- [ ] **[L-2] Make `make test` depend on `make install`** in `Makefile:1-7`.
- [ ] **[L-3] Delete `print(os.environ.get("PATH"))` lines** at module load in `test/viralunity_consensus_test.py:17`, `test/viralunity_meta_test.py:16`, `test/viralunity_meta_cli_test.py:6`. Drop the now-unused `import os` where applicable.

Deferred from quick wins (need design): `CONTRIBUTING.md` content (depends on contribution policy); `[C-2]` proper fix (move pytest/click to dev extras — needs `pyproject.toml`).

### Quick-win verification
- `make test` passes when run after `make install`.
- `viralunity --version` reports `1.1.0` and matches the Dockerfile `LABEL`.
- `grep -rn "AdaptersNotFoundError" viralunity/validators.py` shows it in the import block.
- An Illumina run with a missing `--adapters` value now raises `AdaptersNotFoundError` (the proper exception), not `NameError`. Verify by adding/updating a test that does not mock the validator for this branch — or, at minimum, that `viralunity.validators.AdaptersNotFoundError` is accessible at import time.
- `grep -rn "print(" viralunity/*.py` returns no library hits.
- `LICENSE` reads `2021-2026`.
- `README.md` has installation + usage information visible from the GitHub front page.
- `CHANGELOG.md` exists at repo root.

---

## Phase 2 — High priority backlog

Substantial but well-scoped. One PR per item, in roughly the order below.

- [ ] **[C-2] / [H-3] Split runtime vs dev dependencies.** Move `black`, `pytest`, and (eventually) `ruff` out of `setup.py:install_requires` into a `dev` extras section. Easiest path: introduce a `pyproject.toml` with `[project.optional-dependencies] dev = [...]` and reference it from CI / `CONTRIBUTING.md`.
- [ ] **[H-2] Harden CI**:
  - Add `push: branches: [main]` trigger.
  - Add a `lint` job: `black --check viralunity/ test/` and (when configured) `ruff check`.
  - Add a `dryrun` job: run `test/dryrun_configs/create_dryrun_placeholders.sh`, then `pytest test/viralunity_dryrun_test.py` (the dryrun tests use pytest, not unittest).
  - Once `[C-2]` is done, drop the redundant `pip install -e .` step ([M-11]).
- [ ] **[H-5] Add a dedicated nanopore QC rule.** Today `metagenomics_nanopore.smk` inherits paired-end Illumina QC; verify with a dry-run on the supplied nanopore fixtures (`test/dryrun_configs/metagenomics_nanopore.yaml`) and add `rules/qc_nanopore.smk` with NanoFilt / fastp single-end equivalents.
- [ ] **[H-6] Sweep remaining shell blocks for `set -euo pipefail`.** Targets identified: `rules/alignment_illumina.smk:25-31`, `rules/qc_illumina.smk:39-64`, `consensus_nanopore.smk:77-82`, `rules/consensus_illumina.smk:40-55`. Confirm no rule depends on a non-zero exit being absorbed.
- [ ] **[H-7] Replace `subprocess.Popen(..., shell=True)`** in `viralunity/scripts/calculate_assembly_stats.py:11,19` with `subprocess.run([...], check=True, capture_output=True)`. Or replace counting with `pysam`.
- [ ] **[H-8] Extract shared `_run` subprocess helper** from `viralunity_get_databases_cli.py:25` and `viralunity_build_deacon_index_cli.py:13` into a `viralunity/_subprocess.py`.
- [ ] **[H-9] Consolidate duplicate parsers in `viralunity_get_databases_cli.py`**: merge `_parse_data_report` / `_parse_genome_data_report` and `_reformat_protein_fasta` / `_reformat_genome_fasta` into single parametric functions.

---

## Phase 3 — Architectural / strategic backlog

Each requires a small design discussion.

- [ ] **[H-4] Collapse 4-way consensus workflow duplication.** Each pair (`consensus_illumina.smk` ↔ `consensus_illumina_segmented.smk`; ditto nanopore) can be merged once the top-level rules are parametrized over `SEGMENT_WILDCARD` the way the included `rules/*.smk` already are. Preserve output paths to avoid breaking downstream consumers.
- [ ] **[M-3] Rename custom `FileNotFoundError`** to something like `MissingPathError`. Update all importers in `viralunity/validators.py`, `viralunity/viralunity_create_samplesheet.py`, etc.
- [ ] **[M-5] Docstring sweep on Click handlers.** All five CLI modules: add a one-line summary docstring on each `@click.command()` handler so IDE tooling sees them.
- [ ] **[M-6] Add type hints on Click handler signatures.**
- [ ] **[M-9] Externalize hardcoded tool flags** into config keys: `minimap2` flags in `consensus_*.smk:111` (and segmented variants), DIAMOND `--max-target-seqs` in `metagenomics_diamond_reads_illumina.smk:31`, kraken2 `--minimum-hit-group` in `metagenomics_kraken2_reads_illumina.smk:13`.
- [ ] **[L-4] Align dependency pinning** between `environment.yml` and `setup.py`. Pick patch-pin or minor-pin consistently.
- [ ] **[L-5] Add `pyproject.toml`** with `[tool.black]`, `[tool.ruff]`, `[tool.mypy]`, dev extras, modern PEP 621 metadata. Migrate from `setup.py`.
- [ ] **[L-6] Decide on `forceall`/`lock` semantics** in `viralunity_meta.py:129-131` — almost certainly remove `forceall=True` (defeats caching).
- [ ] **Orchestration deduplication** between `viralunity_consensus.py` and `viralunity_meta.py` (`validate_args`/`generate_config_file`/`main` skeletons) — extract a `PipelineRunner` base or functional template.

---

## Phase 4 — Polish

- [ ] **[L-1 cont.]** Add `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `.github/ISSUE_TEMPLATE/`, `.github/PULL_REQUEST_TEMPLATE.md`.
- [ ] **[M-2] Simplify or document** the `args.update(...)` line in `viralunity/viralunity_meta_cli.py:363`.
- [ ] **[M-8] Helper script polish**: rename `input` parameter; add module + function docstrings; add types in `calculate_assembly_stats.py` and `rename_sequences.py`.
- [ ] **[M-12] Tag `v1.1.0`** retroactively; document a release process.
- [ ] **[L-7] Update Docker `LABEL description`** to use `viralunity._description`.

---

## Cross-cutting recommendation

The branch is markedly more mature than `main`: per-rule conda envs, resources declarations, shell safety, dryrun integration tests, sphinx docs, and a clean Click-based CLI surface are all real wins. After Phase 1 lands, the highest-leverage Phase-2 item is **[H-5]** (nanopore QC rule) because the metagenomics-nanopore path likely doesn't work today — but it's hidden by the dry-run tests only verifying DAG structure rather than data correctness. The other Phase-2 items mostly raise the quality floor on already-functional code.

The Critical findings are all small, mechanical fixes — most of the engagement value is in catching them before this branch merges to `main`.
