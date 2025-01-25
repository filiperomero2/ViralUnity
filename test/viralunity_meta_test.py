import unittest
from unittest.mock import patch
from viralunity.viralunity_meta import get_args, validate_args
import os

print(os.environ.get("PATH"))


class Test_RequiredArgs(unittest.TestCase):
    all_args = [
        "--data-type",
        "illumina",
        "--sample-sheet",
        "sample_sheet.csv",
        "--config-file",
        "config_file.yaml",
        "--output",
        "output_dir",
        "--run-name",
        "run_name",
        "--kraken2-database",
        "kraken2_db",
        "--krona-database",
        "krona_db",
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

    def test_required_args_fail_when_only_optional_set(self):
        args = [
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
                get_args(args)

    def test_required_args_success_when_only_required_set(self):
        args = [
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
        readArgs = get_args(args)
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
        readArgs = get_args(args)
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


class Test_ValidateArgs(unittest.TestCase):
    def setUp(self):
        self.args = {
            "data_type": "illumina",
            "sample_sheet": "sample_sheet.csv",
            "config_file": "config_file.yaml",
            "output": "output_dir",
            "run_name": "run_name",
            "kraken2_database": "kraken2_db",
            "krona_database": "krona_db",
            "remove_human_reads": False,
            "remove_unclassified_reads": False,
            "adapters": "adapters.fasta",
            "minimum_read_length": 50,
            "trim": 0,
            "create_config_only": False,
            "threads": 1,
            "threads_total": 1,
        }

    @patch("os.path.isfile", side_effect=[True, False, True])
    @patch("os.path.isdir", side_effect=[False, True, True])
    @patch("viralunity.viralunity_meta.validate_sample_sheet", return_value={})
    def test_validate_args_success(
        self, mock_validate_sample_sheet, mock_isdir, mock_isfile
    ):
        samples = validate_args(self.args)
        self.assertEqual(samples, {})
        mock_validate_sample_sheet.assert_called_once_with(
            "sample_sheet.csv", self.args
        )

    @patch("os.path.isfile", side_effect=[False, True, True, True, True])
    def test_validate_args_sample_sheet_not_exist(self, mock_isfile):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False])
    @patch("os.path.isdir", side_effect=[False, False])
    @patch("viralunity.viralunity_meta.validate_sample_sheet", return_value={})
    def test_validate_args_kraken2_db_not_exist(self, mock_validate_sample_sheet, mock_isdir, mock_isfile):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False])
    @patch("os.path.isdir", side_effect=[False, True, False])
    @patch("viralunity.viralunity_meta.validate_sample_sheet", return_value={})
    def test_validate_args_krona_db_not_exist(self, mock_validate_sample_sheet, mock_isdir, mock_isfile):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False, False])
    @patch("os.path.isdir", side_effect=[False, True, True])
    @patch("viralunity.viralunity_meta.validate_sample_sheet", return_value={})
    def test_validate_args_adapters_not_exist(self, mock_validate_sample_sheet, mock_isdir, mock_isfile):
        with self.assertRaises(SystemExit):
            validate_args(self.args)


if __name__ == "__main__":
    unittest.main()
