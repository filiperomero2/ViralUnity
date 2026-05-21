## Summary

<!-- 1–3 bullets describing what changed and why. -->

-
-

## Related issue

<!-- Closes #<issue-number>, or "n/a" for housekeeping PRs. -->

## Test plan

<!-- Bulleted checklist of things you ran / verified locally. -->

- [ ] `make test` passes
- [ ] `pytest test/viralunity_dryrun_test.py -v` passes
- [ ] `black --check viralunity/ test/` clean
- [ ] `ruff check viralunity/ test/` clean
- [ ] (if user-facing) `viralunity --help` / subcommand `--help` renders correctly
- [ ] (if Snakefile change) DAG verified via `snakemake -n` against the relevant `test/dryrun_configs/*.yaml`

## Notes for reviewers

<!-- Anything that warrants extra attention: subtle behaviour changes,
backwards-compat considerations, follow-up items intentionally
deferred. -->
