import argparse
import unittest
import os
from viralunity.viralunity_meta_cli import fill_arg_parser_meta

print(os.environ.get("PATH"))


class Test_FillArgParserMeta(unittest.TestCase):
    def test_required_args_fail_when_only_optional_set(self):
        args = [
            "meta",
            "--run-name",
            "run_name",
            "--remove-human-reads",
            "--remove-unclassified-reads",
            "--adapters",
            "adapters.fasta",
            "--minimum-read-length",
            "50",
            "--trim",
            "0",
            "--create-config-only",
            "--threads",
            "1",
        ]
        with self.assertRaises(SystemExit):
            parser = argparse.ArgumentParser()
            subparsers = parser.add_subparsers()
            fill_arg_parser_meta(subparsers)
            parser.parse_args(args)

    def test_required_args_success_when_only_required_set(self):
        args = [
            "meta",
            "--data-type",
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
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_meta(subparsers)
        readArgs = vars(parser.parse_args(args))
        self.assertDictContainsSubset(
            {
                "data_type": "illumina",
                "sample_sheet": "sample_sheet.csv",
                "config_file": "config_file.yaml",
                "output": "output_dir",
                "kraken2_database": "kraken2_db",
                "krona_database": "krona_db",
            },
            readArgs,
        )

    def test_default_values_optional_args(self):
        args = [
            "meta",
            "--data-type",
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
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_meta(subparsers)
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
                "kraken2_database": "kraken2_db",
                "krona_database": "krona_db",
                "minimum_read_length": 50,
                "output": "output_dir",
                "remove_human_reads": False,
                "remove_unclassified_reads": False,
                "run_name": "undefined",
                "sample_sheet": "sample_sheet.csv",
                "threads": 1,
                "threads_total": 1,
                "trim": 0,
            },
        )

if __name__ == "__main__":
    unittest.main()
