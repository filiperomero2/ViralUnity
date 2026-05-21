# ViralUnity — Code Review Backlog

**Branch:** `feature/update-meta` (HEAD `0cea653`)
Prioritized checklist of remediation items distilled from `CODE_REVIEW.md`. Severity first, then effort (small → large) within each severity.

**Status:** all four phases landed. The branch is release-ready; tagging `v1.1.0` is the only outstanding maintainer action — see [`RELEASING.md`](RELEASING.md). Per-item commit hashes are recorded inline so each item is traceable.

---

## Phase 1 — Quick wins ✅ DONE

Landed in commit `84a17c0` ("code review quick wins"). Mechanical, low-risk, no design discussion needed.

- [x] **[C-1] Add `AdaptersNotFoundError` to the import list** in `viralunity/validators.py:7-19`. Add `from e` chaining at line 122-125 while we're there. → `84a17c0`
- [x] **[C-3] Bump Dockerfile version label** at `Dockerfile:2` from `1.0.3` to `1.1.0`. → `84a17c0`
- [x] **[M-1] Replace `print(...)` in library code** at `viralunity/viralunity_meta.py:226` with `logger.warning`. → `84a17c0`
- [x] **[M-4] Add `from e` exception chaining** at all `raise X(...)` inside `except ... as e` blocks in `viralunity/validators.py` and `viralunity/config_generator.py`. → `84a17c0`
- [x] **[M-7] Bump `LICENSE`** copyright to `2021-2026`. → `84a17c0`
- [x] **[M-10] Fix deprecated pandas idiom** in `viralunity/scripts/calculate_assembly_stats.py`: replace `delim_whitespace=True` with `sep=r"\s+"` (and resolve the `# check this` comment by removing it). → `84a17c0`
- [x] **[H-1] Restore a minimal usable README**: keep the RTD link, but add 30-50 lines covering install, two-line usage for each subcommand, link to RTD for full docs. → `84a17c0`
- [x] **[L-1] Add minimal `CHANGELOG.md`** documenting the quick wins. → `84a17c0` (expanded in `0cea653`)
- [x] **[L-2] Make `make test` depend on `make install`** in `Makefile:1-7`. → `84a17c0`
- [x] **[L-3] Delete `print(os.environ.get("PATH"))` lines** at module load in `test/viralunity_consensus_test.py:17`, `test/viralunity_meta_test.py:16`, `test/viralunity_meta_cli_test.py:6`. Drop the now-unused `import os` where applicable. → `84a17c0`

---

## Phase 2 — High priority backlog ✅ DONE (H-5 deferred by design)

- [x] **[C-2] / [H-3] Split runtime vs dev dependencies.** Migrated `setup.py` → `pyproject.toml` (PEP 621) with `[project.optional-dependencies] dev = ["black>=24.4", "pytest>=9.0", "ruff>=0.5"]`. → `ea21537`
- [x] **[H-2] Harden CI** — `push: branches: [main]` trigger; dedicated `lint` job (`black --check`, `ruff check`); dedicated `dryrun` job; redundant `pip install -e .` step removed once dev extras worked. → `c59d841`
- [ ] **[H-5] Add a dedicated nanopore QC rule.** **Deferred by design.** Nanopore QC is performed at basecalling time, so a workflow-level QC rule is redundant. Decision recorded in the Phase 2 plan.
- [x] **[H-6] Sweep remaining shell blocks for `set -euo pipefail`.** → `7021e95`
- [x] **[H-7] Replace `subprocess.Popen(..., shell=True)`** in `calculate_assembly_stats.py` with `subprocess.run([...])`. → `e9b239b`
- [x] **[H-8] Extract shared `_run` subprocess helper** into `viralunity/_subprocess.py` (as `run_command`). → `842b82c`
- [x] **[H-9] Consolidate duplicate parsers in `viralunity_get_databases_cli.py`** — `_parse_data_report` / `_parse_genome_data_report` merged into a single parametric `_parse_data_report_jsonl(report_path, key_extractor=...)`. → `6ca373b`

---

## Phase 3 — Architectural / strategic backlog ✅ DONE (H-4 partial)

- [~] **[H-4] Collapse 4-way consensus workflow duplication** — **partial.** Extracted shared rules into `rules/consensus_illumina_common.smk` and `rules/consensus_nanopore_common.smk`. The full 4 → 2 file collapse stays deferred because the segment-specific shell logic in `organize_files` / `summarize_isnvs` / `sanitize_reference` doesn't collapse cleanly. → `c202e94`
- [x] **[M-3] Rename custom `FileNotFoundError`** to `ViralUnityFileNotFoundError`. → `8bfd7c3`
- [x] **[M-5] Docstring sweep on Click handlers.** → `4f46fb8`
- [x] **[M-6] Add type hints on Click handler signatures.** All 18 handlers across the five CLI modules. → `480b7c4`
- [x] **[M-9] Externalize hardcoded tool flags** into config keys: `minimap2_consensus_align_flags`, `diamond_max_target_seqs`, `kraken2_extra_flags`. Defaults preserve current behaviour exactly. → `680049b`
- [x] **[L-4] Align dependency pinning** between `environment.yml` and `pyproject.toml` on `>=X.Y` semantics, with `<8` upper bound on snakemake. → `de5d545`
- [x] **[L-5] Add `pyproject.toml`** with `[tool.black]`, `[tool.ruff]`, `[tool.mypy]`, `[tool.pytest.ini_options]`, dev extras, PEP 621 metadata. (Pre-commit config still deferred.) → `ea21537` + `4dcbd0b`
- [x] **[L-6] Decide on `forceall`/`lock` semantics** — dropped `forceall=True, lock=False, workdir=os.getcwd()` from the meta pipeline's `snakemake(...)` call. → `4f650ec` (and test re-framed in `4151573`)
- [x] **Orchestration deduplication** between `viralunity_consensus.py` and `viralunity_meta.py` — extracted shared `start_config`, `run_workflow`, `run_pipeline` into `viralunity/_orchestrator.py`. Each pipeline still exposes `validate_args` / `generate_config_file` / `run_snakemake_workflow` / `main` so existing test patches keep working. → `1c3d127`

---

## Phase 4 — Polish ✅ DONE (tag pending maintainer)

- [x] **[L-1 cont.]** Added `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` (links to Contributor Covenant 2.1), `.github/ISSUE_TEMPLATE/{bug_report,feature_request}.md`, `.github/PULL_REQUEST_TEMPLATE.md`. → `4eabbd7`
- [x] **[M-2] Simplify or document** the `args.update(...)` lines in `viralunity/viralunity_meta_cli.py` — deleted as proven-dead code (`_build_meta_args` returns the same dict it mutated, so the comprehension is empty by construction). → `ff92ef1`
- [x] **[M-8] Helper script polish** for `calculate_assembly_stats.py` and `rename_sequences.py`: module + function docstrings, type hints, renamed the `input` parameter that shadowed the Python builtin. → `2ec2a4a`
- [~] **[M-12] Tag `v1.1.0`** — **documentation done; maintainer action pending.** CHANGELOG backfilled with all Phase-2/3/4 entries under `[1.1.0] - 2026-05-21`; `RELEASING.md` documents the full release workflow (bump `__version__`, bump Dockerfile LABEL, commit, tag, push). Tag itself stays with the maintainer per Phase 4 plan: `git tag -a v1.1.0 -m "Release v1.1.0" && git push origin v1.1.0`. → `0cea653`
- [x] **[L-7] Update Docker `LABEL description`** to match `viralunity._description` ("A pipeline for viral metagenomics analysis."). → `71d5e36`

---

## Outstanding follow-ups (not in scope of the four phases)

These were flagged inline during execution but deliberately deferred:

- **`pre-commit` config** — was mentioned alongside `[L-5]`. Skipped because the `[tool.black]` and `[tool.ruff]` configs in `pyproject.toml` give the same enforcement when contributors choose to run the tools, and CI gates the same checks on PRs.
- **Full H-4 collapse (4 → 2 top-level consensus snakefiles)** — extracted the cleanly-parametric rules to `rules/consensus_{illumina,nanopore}_common.smk` in `c202e94`. Going further would require Python-driven shell-script selection that just moves the complexity around.
- **`v1.1.0` git tag** — explicitly left to the maintainer per the Phase-4 plan. `RELEASING.md` has the exact commands.
- **Local-env `pulp<3` workaround** — Snakemake 7.32 uses `pulp.list_solvers` which was renamed to `pulp.listSolvers` in pulp 3.x. Anyone whose env pulls `pulp>=3` needs `pip install 'pulp<3'` for the dryrun tests to run. Not project-pinned because `environment.yml` doesn't lock pulp explicitly; CI uses the conda env which currently resolves to a compatible version.

---

## Cross-cutting recommendation (historical, kept for context)

The branch is markedly more mature than `main`: per-rule conda envs, resources declarations, shell safety, dryrun integration tests, sphinx docs, and a clean Click-based CLI surface are all real wins. After Phase 1 lands, the highest-leverage Phase-2 item is **[H-5]** (nanopore QC rule) because the metagenomics-nanopore path likely doesn't work today — but it's hidden by the dry-run tests only verifying DAG structure rather than data correctness. The other Phase-2 items mostly raise the quality floor on already-functional code.

> **2026-05-21 update:** `[H-5]` was deferred during Phase 2 — nanopore QC happens at basecalling time, so a workflow-level rule is redundant. Everything else listed above has landed.
