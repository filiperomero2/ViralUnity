import argparse
import unittest
from viralunity.viralunity_consensus_cli import fill_arg_parser_consensus

class Test_FillArgParserConsensus(unittest.TestCase):

    def test_get_args_required(self):
        args = [
            "consensus",
            "--run-name",
            "test_run",
            "--primer-scheme",
            "scheme",
            "--minimum-coverage",
            "20",
            "--adapters",
            "--minimum-read-length",
            "50",
            "--trim",
            "0",
            "--create-config-only",
            "--threads",
            "1",
            "--threads-total",
            "1",
        ]
        with self.assertRaises(SystemExit):
            parser = argparse.ArgumentParser()
            subparsers = parser.add_subparsers()
            fill_arg_parser_consensus(subparsers)
            parser.parse_args(args)

    def test_required_args_success_with_reference(self):
        args = [
            "consensus",
            "--data-type",
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
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_consensus(subparsers)
        readArgs = vars(parser.parse_args(args))
        expected_subset = {
            "data_type": "illumina",
            "sample_sheet": "sample_sheet.csv",
            "config_file": "config_file.yaml",
            "output": "output_dir",
            "reference": "reference.fasta",
            "segmented_reference": None,
        }
        for key, value in expected_subset.items():
            self.assertEqual(readArgs[key], value, f"Mismatch for key '{key}'")

    def test_required_args_success_with_segmented_reference(self):
        args = [
            "consensus",
            "--data-type",
            "nanopore",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--segmented-reference",
            "S=/path/to/S.fasta",
            "L=/path/to/L.fasta",
            "M=/path/to/M.fasta",
        ]
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_consensus(subparsers)
        readArgs = vars(parser.parse_args(args))
        self.assertEqual(readArgs["reference"], None)
        self.assertEqual(
            readArgs["segmented_reference"],
            ["S=/path/to/S.fasta", "L=/path/to/L.fasta", "M=/path/to/M.fasta"],
        )

    def test_required_args_success_without_any_reference(self):
        """Both --reference and --segmented-reference are optional at argparse level.
        Mutual exclusion is enforced at validation time."""
        args = [
            "consensus",
            "--data-type",
            "illumina",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
        ]
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_consensus(subparsers)
        readArgs = vars(parser.parse_args(args))
        self.assertIsNone(readArgs["reference"])
        self.assertIsNone(readArgs["segmented_reference"])

    def test_default_values_optional_args(self):
        args = [
            "consensus",
            "--data-type",
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
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_consensus(subparsers)
        readArgs = vars(parser.parse_args(args))
        
        self.assertIn("func", readArgs)
        self.assertIsNotNone(readArgs["func"])
        readArgs.pop("func")  # Remove 'func' before comparison
        
        self.assertDictEqual(
            readArgs,
            {
                "adapters": None,
                "config_file": "config_file.yaml",
                "create_config_only": False,
                "data_type": "illumina",
                "minimum_coverage": 20,
                "minimum_read_length": 50,
                "output": "output_dir",
                "primer_scheme": None,
                "reference": "reference.fasta",
                "segmented_reference": None,
                "run_name": "undefined",
                "sample_sheet": "sample_sheet.csv",
                "threads": 1,
                "threads_total": 1,
                "trim_head": 0,
                "trim_tail": 0,
                "cut_front_mean_quality": 10,
                "cut_tail_mean_quality": 10,
                "cut_right_window_size": 4,
                "cut_right_mean_quality": 15,
                "af_threshold": 0.51,
                "af_isnv_threshold": 0,
                "run_isnv": False,
                "chunk_size": 10000,
                "clair3_model": "r1041_e82_400bps_sup_v500",
                "variant_quality": 20,
                "variant_depth": 10,
                "minimum_map_quality": 30,
            },
        )


if __name__ == "__main__":
    unittest.main()