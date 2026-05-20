# Changelog

All notable changes to this project will be documented in this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aspires to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] — feature/update-meta

### Fixed
- **Critical**: `viralunity/validators.py` raised `AdaptersNotFoundError`
  without importing it; any Illumina run providing a non-existent `--adapters`
  path triggered a `NameError` instead of the intended `AdaptersNotFoundError`.
  Existing tests mocked the validator so this did not surface in CI.
- `Dockerfile` `LABEL version` updated from `1.0.3` to `1.1.0` to match
  `viralunity/__init__.py:__version__`.
- Replaced `print("No samples were provided.")` in
  `viralunity/viralunity_meta.py` with `logger.warning(...)`.
- Added `from e` exception chaining to all custom-exception re-raises inside
  `except ... as e` blocks in `viralunity/validators.py` and
  `viralunity/config_generator.py`.
- Replaced deprecated `pd.read_csv(..., delim_whitespace=True)` with
  `sep=r"\s+"` in `viralunity/scripts/python/calculate_assembly_stats.py`.
- `make test` now depends on `make install` so a fresh clone runs cleanly.
- Removed leftover `print(os.environ.get("PATH"))` debug calls from
  `test/viralunity_consensus_test.py` and `test/viralunity_meta_test.py`.

### Changed
- LICENSE copyright extended to `2021-2026`.
- README restored to a usable shape: install + quick-start blocks for each
  subcommand, with the ReadTheDocs link kept as the canonical source.

## [1.1.0]
See git history for the substantial changes on this branch versus `main`:
per-rule conda environments, resources declarations, set -euo pipefail
sweep, dryrun integration tests, sphinx docs, argparse→Click CLI migration,
new `get-databases` / `build-deacon-index` subcommands, metagenomics
expansion (DIAMOND, reference assembly, dehosting, kraken2 reads/contigs).
