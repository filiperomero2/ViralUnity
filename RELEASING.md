# Releasing ViralUnity

This document describes how to cut a versioned release. Only project
maintainers need to read it.

## Versioning

ViralUnity follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
The version string is single-sourced from
[`viralunity/__init__.py`](viralunity/__init__.py): the value of
`__version__` is read dynamically by `pyproject.toml` through
`[tool.setuptools.dynamic]`, so updating it in one place propagates to
the Python package metadata.

The `Dockerfile`'s `LABEL version="X.Y.Z"` is not derived
automatically — it must be bumped in lockstep.

## Cutting a release (X.Y.Z)

1. **Bump `__version__`** in `viralunity/__init__.py` to `X.Y.Z`.

2. **Bump the Dockerfile LABEL** at `Dockerfile:2` to
   `LABEL version="X.Y.Z"`. Keep the `description` line aligned with
   `viralunity._description` (this should already match).

3. **Update `CHANGELOG.md`**: rename the most recent unreleased section
   (or add a new one) to `## [X.Y.Z] - YYYY-MM-DD`, grouping entries
   under `### Added` / `### Changed` / `### Fixed` / `### Refactored`.

4. **Commit** the three edits as a single commit:

   ```bash
   git add viralunity/__init__.py Dockerfile CHANGELOG.md
   git commit -m "release: X.Y.Z"
   ```

5. **Tag**:

   ```bash
   git tag -a vX.Y.Z -m "Release vX.Y.Z"
   ```

6. **Push** the branch and the tag:

   ```bash
   git push origin main
   git push origin vX.Y.Z
   ```

7. **Build and publish the Docker image** (if applicable):

   ```bash
   docker build -t institutotodospelasaude/viralunity:X.Y.Z .
   docker tag institutotodospelasaude/viralunity:X.Y.Z \
              institutotodospelasaude/viralunity:latest
   docker push institutotodospelasaude/viralunity:X.Y.Z
   docker push institutotodospelasaude/viralunity:latest
   ```

## Verifying

After the tag is pushed:

```bash
# Python package version
viralunity --version
# -> viralunity, version X.Y.Z

# Docker image labels
docker inspect institutotodospelasaude/viralunity:X.Y.Z \
  | jq '.[0].Config.Labels'
# -> { "version": "X.Y.Z", "description": "A pipeline for viral metagenomics analysis." }
```

`viralunity.__version__`, the Dockerfile `LABEL version`, the git tag,
and the topmost `CHANGELOG.md` section should all agree.

## Hotfix releases

For a patch release on top of a tag (e.g. `v1.1.0` → `v1.1.1`):

1. Branch from the tag: `git checkout -b release/1.1.1 v1.1.0`.
2. Apply the fix, then follow steps 1–7 above with `X.Y.Z = 1.1.1`.
3. Open a PR back into `main` so the fix lands there too.

## Backfilling a missed tag

If a release commit landed on `main` but no tag was created, the tag
can be added retroactively at that commit:

```bash
git tag -a vX.Y.Z <commit-sha> -m "Release vX.Y.Z"
git push origin vX.Y.Z
```
