import unittest
from unittest.mock import patch
from viralunity.viralunity_meta import (
    validate_args,
    generate_config_file,
    main,
)
import os

print(os.environ.get("PATH"))

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
    def test_validate_args_kraken2_db_not_exist(
        self, mock_validate_sample_sheet, mock_isdir, mock_isfile
    ):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False])
    @patch("os.path.isdir", side_effect=[False, True, False])
    @patch("viralunity.viralunity_meta.validate_sample_sheet", return_value={})
    def test_validate_args_krona_db_not_exist(
        self, mock_validate_sample_sheet, mock_isdir, mock_isfile
    ):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False, False])
    @patch("os.path.isdir", side_effect=[False, True, True])
    @patch("viralunity.viralunity_meta.validate_sample_sheet", return_value={})
    def test_validate_args_adapters_not_exist(
        self, mock_validate_sample_sheet, mock_isdir, mock_isfile
    ):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, True])
    @patch("viralunity.viralunity_meta.validate_sample_sheet", return_value={})
    def test_validate_args_config_file_exists(
        self, mock_validate_sample_sheet, mock_isfile
    ):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False])
    @patch("os.path.isdir", side_effect=[True])
    @patch("viralunity.viralunity_meta.validate_sample_sheet", return_value={})
    def test_validate_args_output_dir_exists(
        self, mock_validate_sample_sheet, mock_isdir, mock_isfile
    ):
        with self.assertRaises(SystemExit):
            validate_args(self.args)


class Test_GenerateConfigFile(unittest.TestCase):
    def setUp(self):
        self.args = {
            "sample_sheet": "sample_sheet.csv",
            "config_file": "config_file.yaml",
            "output": "output_dir",
            "run_name": "run_name",
            "kraken2_database": "kraken2_db",
            "krona_database": "krona_db",
            "remove_human_reads": True,
            "remove_unclassified_reads": False,
            "adapters": "adapters.fasta",
            "minimum_read_length": 50,
            "trim": 0,
            "create_config_only": False,
            "threads": 1,
            "threads_total": 1,
        }

    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    @patch("viralunity.viralunity_meta.define_job_id", return_value="job_id")
    def test_generate_config_file_illumina(self, mock_define_job_id, mock_open):
        self.args["data_type"] = "illumina"
        self.samples = {
            "sample1": ["sample1_R1.fastq", "sample1_R2.fastq"],
            "sample2": ["sample2_R1.fastq", "sample2_R2.fastq"],
        }
        generate_config_file(self.samples, self.args)
        mock_open.assert_called_once_with("config_file.yaml", "w")
        handle = mock_open()
        handle.write.assert_any_call("samples:\n")
        handle.write.assert_any_call(
            "    sample-sample1: sample1_R1.fastq sample1_R2.fastq\n"
        )
        handle.write.assert_any_call(
            "    sample-sample2: sample2_R1.fastq sample2_R2.fastq\n"
        )
        handle.write.assert_any_call("data: illumina\n")
        handle.write.assert_any_call("kraken2_database: kraken2_db\n")
        handle.write.assert_any_call("krona_database: krona_db\n")
        handle.write.assert_any_call("threads: 1\n")
        handle.write.assert_any_call("output: output_dir/job_id/\n")
        handle.write.assert_any_call("adapters: adapters.fasta\n")
        handle.write.assert_any_call("minimum_length: 50\n")
        handle.write.assert_any_call("trim: 0\n")
        handle.write.assert_any_call("remove_human_reads: True\n")
        handle.write.assert_any_call("remove_unclassified_reads: False\n")

    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    @patch("viralunity.viralunity_meta.define_job_id", return_value="job_id")
    def test_generate_config_file_nanopore(self, mock_define_job_id, mock_open):
        self.args["data_type"] = "nanopore"
        self.samples = {
            "sample1": ["sample1.fastq"],
            "sample2": ["sample2.fastq"],
        }
        generate_config_file(self.samples, self.args)
        mock_open.assert_called_once_with("config_file.yaml", "w")
        handle = mock_open()
        handle.write.assert_any_call("samples:\n")
        handle.write.assert_any_call("    sample-sample1: sample1.fastq\n")
        handle.write.assert_any_call("    sample-sample2: sample2.fastq\n")
        handle.write.assert_any_call("data: nanopore\n")
        handle.write.assert_any_call("kraken2_database: kraken2_db\n")
        handle.write.assert_any_call("krona_database: krona_db\n")
        handle.write.assert_any_call("threads: 1\n")
        handle.write.assert_any_call("output: output_dir/job_id/\n")
        handle.write.assert_any_call("remove_human_reads: True\n")
        handle.write.assert_any_call("remove_unclassified_reads: False\n")


class Test_MainFunction(unittest.TestCase):
    @patch("viralunity.viralunity_meta.validate_args", return_value={})
    @patch("viralunity.viralunity_meta.generate_config_file")
    @patch("viralunity.viralunity_meta.snakemake", return_value=True)
    def test_main_success(
        self,
        mock_snakemake,
        mock_generate_config_file,
        mock_validate_args,
    ):
        result = main({
            "config_file": "config_file.yaml",
            "threads_total": 1,
            "data_type": "illumina",
            "create_config_only": False,
        })
        self.assertEqual(result, 0)
        mock_snakemake.assert_called_once()

    @patch("viralunity.viralunity_meta.validate_args", return_value={})
    @patch("viralunity.viralunity_meta.generate_config_file")
    @patch("viralunity.viralunity_meta.snakemake", return_value=True)
    def test_main_create_config_only(
        self,
        mock_snakemake,
        mock_generate_config_file,
        mock_validate_args,
    ):
        result = main({
            "config_file": "config_file.yaml",
            "threads_total": 1,
            "data_type": "illumina",
            "create_config_only": True,
        })
        self.assertEqual(result, 0)
        mock_snakemake.assert_not_called()

    @patch("viralunity.viralunity_meta.validate_args", return_value={})
    @patch("viralunity.viralunity_meta.generate_config_file")
    @patch("viralunity.viralunity_meta.snakemake", return_value=False)
    def test_main_failure(
        self,
        mock_snakemake,
        mock_generate_config_file,
        mock_validate_args,
    ):
        result = main({
            "config_file": "config_file.yaml",
            "threads_total": 1,
            "data_type": "illumina",
            "create_config_only": False,
        })
        self.assertEqual(result, 1)
        mock_snakemake.assert_called_once()


if __name__ == "__main__":
    unittest.main()
