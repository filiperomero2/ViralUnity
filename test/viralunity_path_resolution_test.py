"""Tests for path-argument resolution in the viralunity CLI.

The meta CLI used to pass ``workdir=os.path.dirname(<config_file>)`` to
Snakemake, which silently re-rooted every relative path the user provided on
the command line. The fix is twofold:

1. ``viralunity_meta.main`` and ``viralunity_consensus.main`` now resolve
   every path-typed CLI argument to an absolute path against ``os.getcwd()``
   *before* validation runs.
2. ``viralunity_meta.run_snakemake_workflow`` does not pass ``workdir`` to
   Snakemake, so Snakemake uses the shell's cwd (which is what we want)
   rather than the config file's parent directory.

These tests pin down both behaviours.
"""

import os
import tempfile
import unittest
from unittest.mock import patch

from viralunity import viralunity_consensus, viralunity_meta
from viralunity.validators import (
    META_PATH_ARG_KEYS,
    _is_path_sentinel,
    resolve_path_args,
)

# ---------------------------------------------------------------------------
# _is_path_sentinel
# ---------------------------------------------------------------------------


class TestIsPathSentinel(unittest.TestCase):
    def test_none_is_sentinel(self):
        self.assertTrue(_is_path_sentinel(None))

    def test_empty_string_is_sentinel(self):
        self.assertTrue(_is_path_sentinel(""))
        self.assertTrue(_is_path_sentinel("   "))

    def test_NA_is_sentinel(self):
        self.assertTrue(_is_path_sentinel("NA"))
        self.assertTrue(_is_path_sentinel("  NA  "))

    def test_lowercase_na_is_not_sentinel(self):
        """Matches the case-sensitive ``"NA"`` check used by the validators."""
        self.assertFalse(_is_path_sentinel("na"))

    def test_non_string_is_sentinel(self):
        self.assertTrue(_is_path_sentinel(0))
        self.assertTrue(_is_path_sentinel(False))
        self.assertTrue(_is_path_sentinel(123))

    def test_real_paths_are_not_sentinel(self):
        self.assertFalse(_is_path_sentinel("/abs/path"))
        self.assertFalse(_is_path_sentinel("relative/path"))
        self.assertFalse(_is_path_sentinel("file.fasta"))


# ---------------------------------------------------------------------------
# resolve_path_args
# ---------------------------------------------------------------------------


class TestResolvePathArgs(unittest.TestCase):
    def test_relative_paths_become_absolute_against_base_dir(self):
        args = {
            "kraken2_database": "databases/kraken2/",
            "host_reference": "databases/host/host.fasta",
        }
        resolve_path_args(args, META_PATH_ARG_KEYS, base_dir="/work/proj")
        self.assertEqual(args["kraken2_database"], "/work/proj/databases/kraken2")
        self.assertEqual(args["host_reference"], "/work/proj/databases/host/host.fasta")

    def test_absolute_paths_are_left_untouched(self):
        args = {"host_reference": "/some/abs/path/host.fasta"}
        resolve_path_args(args, META_PATH_ARG_KEYS, base_dir="/work/proj")
        self.assertEqual(args["host_reference"], "/some/abs/path/host.fasta")

    def test_sentinels_are_left_untouched(self):
        args = {
            "host_reference": "NA",
            "deacon_index": "  NA  ",
            "adapters": None,
            "diamond_database": "",
        }
        resolve_path_args(args, META_PATH_ARG_KEYS, base_dir="/work/proj")
        self.assertEqual(args["host_reference"], "NA")
        self.assertEqual(args["deacon_index"], "  NA  ")
        self.assertIsNone(args["adapters"])
        self.assertEqual(args["diamond_database"], "")

    def test_missing_keys_are_ignored(self):
        args = {"kraken2_database": "databases/k2"}
        resolve_path_args(args, META_PATH_ARG_KEYS, base_dir="/work/proj")
        self.assertEqual(args["kraken2_database"], "/work/proj/databases/k2")
        for k in META_PATH_ARG_KEYS:
            if k == "kraken2_database":
                continue
            self.assertNotIn(k, args)

    def test_segmented_reference_dict_values_are_resolved(self):
        args = {
            "reference": {
                "L": "refs/L.fasta",
                "S": "/abs/refs/S.fasta",
            }
        }
        resolve_path_args(args, ("reference",), base_dir="/work/proj")
        self.assertEqual(
            args["reference"],
            {"L": "/work/proj/refs/L.fasta", "S": "/abs/refs/S.fasta"},
        )

    def test_defaults_to_cwd_when_base_dir_omitted(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_real = os.path.realpath(tmp)
            old_cwd = os.getcwd()
            try:
                os.chdir(tmp_real)
                args = {"output": "results/foo"}
                resolve_path_args(args, ("output",))
                self.assertEqual(args["output"], os.path.join(tmp_real, "results/foo"))
            finally:
                os.chdir(old_cwd)

    def test_returns_args_dict_for_chaining(self):
        args = {"output": "results"}
        returned = resolve_path_args(args, ("output",), base_dir="/x")
        self.assertIs(returned, args)


# ---------------------------------------------------------------------------
# Integration: viralunity_meta.main calls resolve_path_args before validation
# ---------------------------------------------------------------------------


class TestMetaMainResolvesPaths(unittest.TestCase):
    """``main()`` must resolve paths *before* validation so validators
    receive absolute paths."""

    def test_main_resolves_paths_before_validate_args(self):
        captured = {}

        def fake_validate(args):
            captured.update(args)
            return {"sample1": ["/abs/sample1.fastq"]}

        with patch("viralunity.viralunity_meta.validate_args", side_effect=fake_validate), patch(
            "viralunity.viralunity_meta.generate_config_file"
        ), patch(
            "viralunity.viralunity_meta.run_snakemake_workflow", return_value=True
        ), tempfile.TemporaryDirectory() as tmp:
            tmp_real = os.path.realpath(tmp)
            old_cwd = os.getcwd()
            os.chdir(tmp_real)
            try:
                exit_code = viralunity_meta.main(
                    {
                        "data_type": "nanopore",
                        "config_file": "scratch/run.yml",
                        "output": "results/run",
                        "kraken2_database": "databases/k2/",
                        "krona_database": "databases/krona/taxonomy",
                        "taxdump": "databases/taxdump/",
                        "host_reference": "databases/host/host.fasta",
                        "diamond_database": "databases/diamond/viral.dmnd",
                        "taxids": "databases/diamond/protein2taxid.tsv",
                        "deacon_index": "NA",
                        "threads_total": 1,
                        "create_config_only": True,
                    }
                )
            finally:
                os.chdir(old_cwd)

        self.assertEqual(exit_code, 0)
        # Every path-typed arg the test supplied is now absolute and
        # rooted at the temporary cwd we created.
        for key in (
            "config_file",
            "output",
            "kraken2_database",
            "krona_database",
            "taxdump",
            "host_reference",
            "diamond_database",
            "taxids",
        ):
            self.assertTrue(
                os.path.isabs(captured[key]),
                f"{key} should be absolute, got {captured[key]!r}",
            )
            self.assertTrue(
                captured[key].startswith(tmp_real),
                f"{key} should be rooted under {tmp_real}, got {captured[key]!r}",
            )
        # Sentinels are preserved untouched.
        self.assertEqual(captured["deacon_index"], "NA")


class TestMetaWorkdirUsesCwd(unittest.TestCase):
    """``run_snakemake_workflow`` must not pin the working directory to the
    config file's parent directory (this was the bug).

    The fix relies on Snakemake's default behaviour: when ``workdir`` is not
    passed to ``snakemake()``, Snakemake does not ``chdir`` and the shell's
    cwd is used as-is. This test asserts that ``workdir`` is either omitted
    or set to cwd — never the config file's directory.
    """

    def test_workdir_is_not_config_dirname(self):
        captured = {}

        def fake_snakemake(*args, **kwargs):
            captured.update(kwargs)
            return True

        with patch(
            "viralunity.viralunity_meta.snakemake", side_effect=fake_snakemake
        ), tempfile.TemporaryDirectory() as tmp:
            tmp_real = os.path.realpath(tmp)
            old_cwd = os.getcwd()
            os.chdir(tmp_real)
            try:
                viralunity_meta.run_snakemake_workflow(
                    {
                        "data_type": "nanopore",
                        # Pretend the user put the config in a sub-directory
                        # so the old code would have set workdir=that subdir.
                        "config_file": os.path.join(tmp_real, "scratch", "run.yml"),
                        "threads_total": 1,
                    }
                )
            finally:
                os.chdir(old_cwd)

        config_dir = os.path.join(tmp_real, "scratch")
        passed_workdir = captured.get("workdir")
        self.assertNotEqual(
            passed_workdir,
            config_dir,
            "workdir must NOT be the config file's parent directory",
        )
        # Either workdir was omitted (Snakemake defaults to cwd) or it was
        # explicitly set to cwd. Both are acceptable.
        if passed_workdir is not None:
            self.assertEqual(passed_workdir, tmp_real)


# ---------------------------------------------------------------------------
# Integration: viralunity_consensus.main calls resolve_path_args
# ---------------------------------------------------------------------------


class TestConsensusMainResolvesPaths(unittest.TestCase):
    def test_main_resolves_paths_before_validate_args(self):
        captured = {}

        def fake_validate(args):
            captured.update(args)
            return {"sample1": ["/abs/sample1.fastq"]}

        with patch(
            "viralunity.viralunity_consensus.validate_args", side_effect=fake_validate
        ), patch("viralunity.viralunity_consensus.generate_config_file"), patch(
            "viralunity.viralunity_consensus.run_snakemake_workflow", return_value=True
        ), tempfile.TemporaryDirectory() as tmp:
            tmp_real = os.path.realpath(tmp)
            old_cwd = os.getcwd()
            os.chdir(tmp_real)
            try:
                exit_code = viralunity_consensus.main(
                    {
                        "data_type": "illumina",
                        "config_file": "scratch/cons.yml",
                        "output": "out/cons",
                        "reference": "refs/ref.fasta",
                        "primer_scheme": "schemes/scheme.bed",
                        "adapters": "adapters/illumina.fa",
                        "create_config_only": True,
                        "threads_total": 1,
                    }
                )
            finally:
                os.chdir(old_cwd)

        self.assertEqual(exit_code, 0)
        for key in (
            "config_file",
            "output",
            "reference",
            "primer_scheme",
            "adapters",
        ):
            self.assertTrue(
                os.path.isabs(captured[key]),
                f"{key} should be absolute, got {captured[key]!r}",
            )
            self.assertTrue(captured[key].startswith(tmp_real))


if __name__ == "__main__":
    unittest.main()
