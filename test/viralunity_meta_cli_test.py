"""Tests for viralunity meta CLI (click-based)."""

import unittest
from unittest.mock import patch

from click.testing import CliRunner

from viralunity.viralunity_meta_cli import meta


class Test_MetaIlluminaCommand(unittest.TestCase):
    """Tests for `viralunity meta illumina`."""

    def setUp(self):
        self.runner = CliRunner()
        self._required = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--kraken2-database",
            "kraken2_db",
            "--krona-database",
            "krona_db",
        ]

    def _invoke(self, extra_args=None):
        args = self._required + (extra_args or [])
        with patch("viralunity.viralunity_meta_cli.meta_main", return_value=0):
            return self.runner.invoke(meta, args, catch_exceptions=False)

    def test_required_args_missing_causes_error(self):
        """Missing required args should exit with non-zero code."""
        result = self.runner.invoke(meta, ["illumina"])
        self.assertNotEqual(result.exit_code, 0)

    def test_required_args_success(self):
        result = self._invoke()
        self.assertEqual(result.exit_code, 0, result.output)

    def test_default_values_optional_args(self):
        """Check that all optional Illumina meta args have correct defaults."""
        with patch("viralunity.viralunity_meta_cli.meta_main", return_value=0) as mock_main:
            result = self.runner.invoke(meta, self._required, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        args = mock_main.call_args[0][0]
        self.assertEqual(args["data_type"], "illumina")
        self.assertEqual(args["run_name"], "undefined")
        self.assertEqual(args["adapters"], "NA")
        self.assertEqual(args["minimum_read_length"], 50)
        self.assertIsNone(args["trim_head"])
        self.assertIsNone(args["trim_tail"])
        self.assertFalse(args["remove_human_reads"])
        self.assertFalse(args["remove_unclassified_reads"])
        self.assertEqual(args["threads"], 1)
        self.assertEqual(args["threads_total"], 1)
        self.assertFalse(args["create_config_only"])

    def test_remove_human_reads_flag(self):
        with patch("viralunity.viralunity_meta_cli.meta_main", return_value=0) as mock_main:
            result = self.runner.invoke(
                meta, self._required + ["--remove-human-reads"], catch_exceptions=False
            )
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertTrue(mock_main.call_args[0][0]["remove_human_reads"])


class Test_MetaNanoporeCommand(unittest.TestCase):
    """Tests for `viralunity meta nanopore`."""

    def setUp(self):
        self.runner = CliRunner()
        self._required = [
            "nanopore",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--kraken2-database",
            "kraken2_db",
            "--krona-database",
            "krona_db",
        ]

    def test_required_args_success(self):
        with patch("viralunity.viralunity_meta_cli.meta_main", return_value=0):
            result = self.runner.invoke(meta, self._required, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)

    def test_default_values_optional_args(self):
        """Check nanopore-specific defaults for the meta command."""
        with patch("viralunity.viralunity_meta_cli.meta_main", return_value=0) as mock_main:
            result = self.runner.invoke(meta, self._required, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        args = mock_main.call_args[0][0]
        self.assertEqual(args["data_type"], "nanopore")
        self.assertFalse(args["run_polish_racon"])
        self.assertFalse(args["run_polish_medaka"])
        self.assertIsNone(args["medaka_model"])
        self.assertIsNone(args["clair3_model"])

    def test_kraken2_flags(self):
        """Test --no-kraken2-reads disables reads classification."""
        with patch("viralunity.viralunity_meta_cli.meta_main", return_value=0) as mock_main:
            result = self.runner.invoke(
                meta, self._required + ["--no-kraken2-reads"], catch_exceptions=False
            )
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertFalse(mock_main.call_args[0][0]["run_kraken2_reads"])


if __name__ == "__main__":
    unittest.main()
