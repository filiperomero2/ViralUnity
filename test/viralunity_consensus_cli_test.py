"""Tests for viralunity consensus CLI (click-based)."""

import unittest
from unittest.mock import patch
from click.testing import CliRunner
from viralunity.viralunity_consensus_cli import consensus


class Test_ConsensusIlluminaCommand(unittest.TestCase):
    """Tests for `viralunity consensus illumina`."""

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
            "--reference",
            "reference.fasta",
        ]

    def _invoke(self, extra_args=None):
        args = self._required + (extra_args or [])
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ):
            return self.runner.invoke(consensus, args, catch_exceptions=False)

    def test_required_args_missing_causes_error(self):
        """Missing required args should exit with non-zero code."""
        result = self.runner.invoke(consensus, ["illumina"])
        self.assertNotEqual(result.exit_code, 0)

    def test_required_args_success_with_reference(self):
        result = self._invoke()
        self.assertEqual(result.exit_code, 0, result.output)

    def test_required_args_success_without_reference(self):
        """Reference flags are optional at parse time; validated by core logic."""
        args = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
        ]
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ):
            result = self.runner.invoke(consensus, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)

    def test_required_args_success_with_segmented_reference(self):
        args = [
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--segmented-reference",
            "S=/path/to/S.fasta",
            "--segmented-reference",
            "L=/path/to/L.fasta",
        ]
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ) as mock_main:
            result = self.runner.invoke(consensus, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        called_args = mock_main.call_args[0][0]
        self.assertEqual(
            called_args["segmented_reference"],
            {
                "S": "/path/to/S.fasta",
                "L": "/path/to/L.fasta",
            },
        )

    def test_default_values_optional_args(self):
        """Check that all optional args have correct defaults for illumina."""
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ) as mock_main:
            result = self.runner.invoke(
                consensus, self._required, catch_exceptions=False
            )
        self.assertEqual(result.exit_code, 0, result.output)
        args = mock_main.call_args[0][0]
        self.assertEqual(args["data_type"], "illumina")
        self.assertEqual(args["run_name"], "undefined")
        self.assertIsNone(args["adapters"])
        self.assertEqual(args["trim_head"], 0)
        self.assertEqual(args["trim_tail"], 0)
        self.assertEqual(args["cut_front_mean_quality"], 10)
        self.assertEqual(args["cut_tail_mean_quality"], 10)
        self.assertEqual(args["cut_right_window_size"], 4)
        self.assertEqual(args["cut_right_mean_quality"], 15)
        self.assertEqual(args["af_threshold"], 0.51)
        self.assertEqual(args["af_isnv_threshold"], 0.0)
        self.assertFalse(args["run_isnv"])
        self.assertEqual(args["minimum_coverage"], 20)
        self.assertEqual(args["minimum_read_length"], 50)
        self.assertEqual(args["threads"], 1)
        self.assertEqual(args["threads_total"], 1)
        self.assertFalse(args["create_config_only"])


class Test_ConsensusNanoporeCommand(unittest.TestCase):
    """Tests for `viralunity consensus nanopore`."""

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
            "--reference",
            "reference.fasta",
        ]

    def test_required_args_success(self):
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ):
            result = self.runner.invoke(
                consensus, self._required, catch_exceptions=False
            )
        self.assertEqual(result.exit_code, 0, result.output)

    def test_default_values_optional_args(self):
        """Check nanopore-specific defaults."""
        with patch(
            "viralunity.viralunity_consensus_cli.consensus_main", return_value=0
        ) as mock_main:
            result = self.runner.invoke(
                consensus, self._required, catch_exceptions=False
            )
        self.assertEqual(result.exit_code, 0, result.output)
        args = mock_main.call_args[0][0]
        self.assertEqual(args["data_type"], "nanopore")
        self.assertEqual(args["af_threshold"], 0.51)
        self.assertEqual(args["chunk_size"], 10000)
        self.assertEqual(args["clair3_model"], "r1041_e82_400bps_sup_v500")
        self.assertEqual(args["variant_quality"], 20)
        self.assertEqual(args["variant_depth"], 10)
        self.assertEqual(args["minimum_map_quality"], 30)


if __name__ == "__main__":
    unittest.main()
