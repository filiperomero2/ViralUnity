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

    def test_required_args_success_when_only_required_set(self):
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
        self.assertDictContainsSubset(
            {
                "data_type": "illumina",
                "sample_sheet": "sample_sheet.csv",
                "config_file": "config_file.yaml",
                "output": "output_dir",
                "reference": "reference.fasta",
            },
            readArgs,
        )

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
                "run_name": "undefined",
                "sample_sheet": "sample_sheet.csv",
                "threads": 1,
                "threads_total": 1,
                "trim": 0,
            },
        )


if __name__ == "__main__":
    unittest.main()