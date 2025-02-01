import unittest
from unittest.mock import patch
import sys
import pandas as pd
from viralunity.viralunity_consensus import (
    get_args,
    validate_args,
    validate_sample_sheet,
    define_job_id,
    generate_config_file,
    main,
)
import os

print(os.environ.get("PATH"))


class Test_RequiredArgs(unittest.TestCase):

    def test_get_args_required(self):
        args = [
            "--rub-name",
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
            "--reference",
            "reference.fasta",
        ]
        readArgs = get_args(args)
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
        readArgs = get_args(args)
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


class Test_ValidateArgs(unittest.TestCase):
    def setUp(self):
        self.args = {
            "data_type": "illumina",
            "sample_sheet": "sample_sheet.csv",
            "config_file": "config_file.yaml",
            "output": "output_dir",
            "run_name": "run_name",
            "reference": "reference.fasta",
            "primer_scheme": "scheme",
            "minimum_coverage": 5,
            "adapters": "adapters.fasta",
            "minimum_read_length": 50,
            "trim": 0,
            "create_config_only": False,
            "threads": 1,
            "threads_total": 1,
        }

    @patch("os.path.isfile", side_effect=[True, False, True, True, True])
    @patch("os.path.isdir", side_effect=[False, True, True])
    @patch("viralunity.viralunity_consensus.validate_sample_sheet", return_value={})
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

    @patch("os.path.isfile", side_effect=[True, True])
    @patch("viralunity.viralunity_consensus.validate_sample_sheet", return_value={})
    def test_validate_args_config_file_exists(
        self, mock_validate_sample_sheet, mock_isfile
    ):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False])
    @patch("os.path.isdir", side_effect=[True])
    @patch("viralunity.viralunity_consensus.validate_sample_sheet", return_value={})
    def test_validate_args_output_dir_exists(
        self, mock_validate_sample_sheet, mock_isdir, mock_isfile
    ):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False, False])
    @patch("os.path.isdir", side_effect=[False])
    @patch("viralunity.viralunity_consensus.validate_sample_sheet", return_value={})
    def test_validate_args_reference_file_not_exist(
        self, mock_validate_sample_sheet, mock_isdir, mock_isfile
    ):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False, True, True, True])
    @patch("os.path.isdir", side_effect=[False, True, True])
    @patch("viralunity.viralunity_consensus.validate_sample_sheet", return_value={})
    def test_validate_args_primer_scheme_not_set(
        self, mock_validate_sample_sheet, mock_isdir, mock_isfile
    ):
        self.args["primer_scheme"] = None
        samples = validate_args(self.args)
        self.assertEqual(self.args["primer_scheme"], "NA")

    @patch("os.path.isfile", side_effect=[True, False, True, True, False])
    @patch("os.path.isdir", side_effect=[False, True, True])
    @patch("viralunity.viralunity_consensus.validate_sample_sheet", return_value={})
    def test_validate_args_illumina_adapters_not_exist(
        self, mock_validate_sample_sheet, mock_isdir, mock_isfile
    ):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False, True, True])
    @patch("os.path.isdir", side_effect=[False, True, True])
    @patch("viralunity.viralunity_consensus.validate_sample_sheet", return_value={})
    def test_validate_args_illumina_adapters_not_set(
        self, mock_validate_sample_sheet, mock_isdir, mock_isfile
    ):
        self.args["adapters"] = None
        with self.assertRaises(SystemExit):
            validate_args(self.args)


class Test_ValidateSampleSheet(unittest.TestCase):
    def setUp(self):
        self.samplesheet_dataframe = pd.DataFrame(
            {
                0: ["sample1", "sample2"],
                1: ["R1_sample1.fastq", "R1_sample2.fastq"],
                2: ["R2_sample1.fastq", "R2_sample2.fastq"],
            }
        )

    @patch("os.path.isfile", return_value=True)
    def test_validate_sample_sheet_illumina(self, mock_isfile):
        with patch("pandas.read_csv", return_value=self.samplesheet_dataframe):
            samples = validate_sample_sheet(
                "sample_sheet.csv", {"data_type": "illumina"}
            )
            self.assertEqual(
                samples,
                {
                    "sample1": ["R1_sample1.fastq", "R2_sample1.fastq"],
                    "sample2": ["R1_sample2.fastq", "R2_sample2.fastq"],
                },
            )

    @patch("os.path.isfile", return_value=False)
    def test_not_validate_sample_sheet_illumina(self, mock_isfile):
        with patch("pandas.read_csv", return_value=self.samplesheet_dataframe):
            with self.assertRaises(SystemExit):
                validate_sample_sheet("sample_sheet.csv", {"data_type": "illumina"})

    @patch("os.path.isfile", return_value=True)
    def test_validate_sample_sheet_nanopore(self, mock_isfile):
        with patch("pandas.read_csv", return_value=self.samplesheet_dataframe):
            samples = validate_sample_sheet(
                "sample_sheet.csv", {"data_type": "nanopore"}
            )
            self.assertEqual(
                samples,
                {
                    "sample1": ["R1_sample1.fastq"],
                    "sample2": ["R1_sample2.fastq"],
                },
            )

    @patch("os.path.isfile", return_value=False)
    def test_not_validate_sample_sheet_nanopore(self, mock_isfile):
        with patch("pandas.read_csv", return_value=self.samplesheet_dataframe):
            with self.assertRaises(SystemExit):
                validate_sample_sheet("sample_sheet.csv", {"data_type": "nanopore"})


class Test_DefineJobId(unittest.TestCase):
    def test_define_job_id(self):
        self.assertRegex(
            define_job_id({"run_name": "viral_run_name"}),
            r"job_\d{4}-\d{2}-\d{2}_viral_run_name",
        )


class Test_GenerateConfigFIle(unittest.TestCase):
    def setUp(self):
        self.args = {
            "sample_sheet": "sample_sheet.csv",
            "config_file": "config_file.yaml",
            "output": "output_dir",
            "run_name": "run_name",
            "reference": "reference.fasta",
            "primer_scheme": "scheme",
            "minimum_coverage": 5,
            "adapters": "adapters.fasta",
            "minimum_read_length": 50,
            "trim": 0,
            "create_config_only": False,
            "threads": 1,
            "threads_total": 1,
        }
    
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    @patch("viralunity.viralunity_consensus.define_job_id", return_value="job_id")
    def test_generate_config_file_illumina(self, mock_define_job_id, mock_open):
        self.args['data_type'] = 'illumina'
        self.samples = {
            "sample1": ["R1_sample1.fastq", "R2_sample1.fastq"],
            "sample2": ["R1_sample2.fastq", "R2_sample2.fastq"],
        }
        
        generate_config_file(self.samples, self.args)
        
        mock_open.assert_called_once_with("config_file.yaml", "w")
        handle = mock_open()
        handle.write.assert_any_call("samples:\n")
        handle.write.assert_any_call("    sample-sample1: R1_sample1.fastq R2_sample1.fastq\n")
        handle.write.assert_any_call("    sample-sample2: R1_sample2.fastq R2_sample2.fastq\n")
        handle.write.assert_any_call("data: illumina\n")
        handle.write.assert_any_call("reference: reference.fasta\n")
        handle.write.assert_any_call("scheme: scheme\n")
        handle.write.assert_any_call("minimum_depth: 5\n")
        handle.write.assert_any_call("threads: 1\n")
        handle.write.assert_any_call(f"workflow_path: {sys.path[0]}\n")
        handle.write.assert_any_call("output: output_dir/job_id/\n")
        handle.write.assert_any_call("adapters: adapters.fasta\n")
        handle.write.assert_any_call("minimum_length: 50\n")
        handle.write.assert_any_call("trim: 0\n")
        
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    @patch("viralunity.viralunity_consensus.define_job_id", return_value="job_id")    
    def test_generate_config_file_nanopore(self, mock_define_job_id, mock_open):
        self.args['data_type'] = 'nanopore'
        self.samples = {
            "sample1": ["R1_sample1.fastq"],
            "sample2": ["R1_sample2.fastq"],
        }
        
        generate_config_file(self.samples, self.args)
        
        mock_open.assert_called_once_with("config_file.yaml", "w")
        handle = mock_open()
        handle.write.assert_any_call("samples:\n")
        handle.write.assert_any_call("    sample-sample1: R1_sample1.fastq\n")
        handle.write.assert_any_call("    sample-sample2: R1_sample2.fastq\n")
        handle.write.assert_any_call("data: nanopore\n")
        handle.write.assert_any_call("reference: reference.fasta\n")
        handle.write.assert_any_call("scheme: scheme\n")
        handle.write.assert_any_call("minimum_depth: 5\n")
        handle.write.assert_any_call("threads: 1\n")
        handle.write.assert_any_call(f"workflow_path: {sys.path[0]}\n")
        handle.write.assert_any_call("output: output_dir/job_id/\n")  
        
class Test_MainFunction(unittest.TestCase):
    @patch("viralunity.viralunity_consensus.get_args",
        return_value={
            "config_file": "config_file.yaml",
            "threads_total": 1,
            "data_type": "illumina",
            "create_config_only": False,
        },)
    @patch("viralunity.viralunity_consensus.validate_args", return_value={})
    @patch("viralunity.viralunity_consensus.generate_config_file")
    @patch("viralunity.viralunity_consensus.snakemake", return_value=True)
    def test_main_success(self, mock_snakemake, mock_generate_config_file, mock_validate_args, mock_get_args):
        result = main()
        self.assertEqual(result, 0)
        mock_snakemake.assert_called_once()
        
        
    @patch("viralunity.viralunity_consensus.get_args",
        return_value={
            "config_file": "config_file.yaml",
            "threads_total": 1,
            "data_type": "illumina",
            "create_config_only": True,
        },)
    @patch("viralunity.viralunity_consensus.validate_args", return_value={})
    @patch("viralunity.viralunity_consensus.generate_config_file")
    @patch("viralunity.viralunity_consensus.snakemake", return_value=False)
    def test_main_create_config_only(self, mock_snakemake, mock_generate_config_file, mock_validate_args, mock_get_args):
        result = main()
        self.assertEqual(result, 0)
        mock_snakemake.assert_not_called()

if __name__ == "__main__":
    unittest.main()
