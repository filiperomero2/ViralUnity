import sys
import unittest
from unittest.mock import mock_open, patch

import pandas as pd

from viralunity.exceptions import (
    AdaptersNotFoundError,
    SampleConfigurationNotFoundError,
    ValidationError,
)
from viralunity.viralunity_consensus import (
    generate_config_file,
    main,
    validate_args,
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
            "segmented_reference": None,
            "primer_scheme": "scheme",
            "minimum_coverage": 5,
            "adapters": "adapters.fasta",
            "minimum_read_length": 50,
            "trim": 0,
            "create_config_only": False,
            "threads": 1,
            "threads_total": 1,
        }

    @patch("viralunity.viralunity_consensus.validate_illumina_requirements")
    @patch("viralunity.viralunity_consensus.validate_consensus_requirements")
    @patch("viralunity.viralunity_consensus.get_samples_from_args")
    def test_validate_args_success(
        self, mock_get_samples, mock_validate_consensus, mock_validate_illumina
    ):
        """Test successful validation with all validators passing."""
        mock_get_samples.return_value = {"sample1": ["file1_R1.fastq", "file1_R2.fastq"]}

        samples = validate_args(self.args)

        self.assertIn("sample1", samples)
        mock_get_samples.assert_called_once_with(self.args)
        mock_validate_consensus.assert_called_once_with(self.args)
        mock_validate_illumina.assert_called_once_with(self.args)

    @patch("viralunity.viralunity_consensus.get_samples_from_args")
    def test_validate_args_sample_sheet_not_exist(self, mock_get_samples):
        """Test validation fails when sample sheet cannot be retrieved."""
        mock_get_samples.side_effect = SampleConfigurationNotFoundError("Sample sheet not found")

        with self.assertRaises(SampleConfigurationNotFoundError):
            validate_args(self.args)

        mock_get_samples.assert_called_once_with(self.args)

    @patch("viralunity.viralunity_consensus.validate_illumina_requirements")
    @patch("viralunity.viralunity_consensus.validate_consensus_requirements")
    @patch("viralunity.viralunity_consensus.get_samples_from_args")
    def test_validate_args_config_file_exists(
        self, mock_get_samples, mock_validate_consensus, mock_validate_illumina
    ):
        """Test validation succeeds even if config file already exists."""
        mock_get_samples.return_value = {"sample1": ["file1_R1.fastq", "file1_R2.fastq"]}

        samples = validate_args(self.args)

        self.assertIn("sample1", samples)

    @patch("viralunity.viralunity_consensus.validate_illumina_requirements")
    @patch("viralunity.viralunity_consensus.validate_consensus_requirements")
    @patch("viralunity.viralunity_consensus.get_samples_from_args")
    def test_validate_args_output_dir_exists(
        self, mock_get_samples, mock_validate_consensus, mock_validate_illumina
    ):
        """Test validation succeeds even if output directory already exists."""
        mock_get_samples.return_value = {"sample1": ["file1_R1.fastq", "file1_R2.fastq"]}

        samples = validate_args(self.args)

        self.assertIn("sample1", samples)

    @patch("viralunity.viralunity_consensus.validate_consensus_requirements")
    @patch("viralunity.viralunity_consensus.get_samples_from_args")
    def test_validate_args_reference_file_not_exist(
        self, mock_get_samples, mock_validate_consensus
    ):
        """Test validation fails when reference file doesn't exist."""
        mock_get_samples.return_value = {"sample1": ["file1.fastq"]}
        mock_validate_consensus.side_effect = ValidationError(
            "Reference sequence file does not exist"
        )

        with self.assertRaises(ValidationError):
            validate_args(self.args)

    @patch("viralunity.viralunity_consensus.validate_illumina_requirements")
    @patch("viralunity.viralunity_consensus.validate_consensus_requirements")
    @patch("viralunity.viralunity_consensus.get_samples_from_args")
    def test_validate_args_both_reference_and_segmented_fails(
        self, mock_get_samples, mock_validate_consensus, mock_validate_illumina
    ):
        """Test validation fails when both --reference and --segmented-reference are provided."""
        mock_get_samples.return_value = {"sample1": ["file1.fastq"]}
        mock_validate_consensus.side_effect = ValidationError(
            "--reference and --segmented-reference are mutually exclusive."
        )
        self.args["segmented_reference"] = ["S=ref_S.fasta"]

        with self.assertRaises(ValidationError):
            validate_args(self.args)

    @patch("viralunity.viralunity_consensus.validate_illumina_requirements")
    @patch("viralunity.viralunity_consensus.validate_consensus_requirements")
    @patch("viralunity.viralunity_consensus.get_samples_from_args")
    def test_validate_args_no_reference_fails(
        self, mock_get_samples, mock_validate_consensus, mock_validate_illumina
    ):
        """Test validation fails when neither --reference nor --segmented-reference are provided."""
        mock_get_samples.return_value = {"sample1": ["file1.fastq"]}
        mock_validate_consensus.side_effect = ValidationError("A reference is required.")
        self.args["reference"] = None
        self.args["segmented_reference"] = None

        with self.assertRaises(ValidationError):
            validate_args(self.args)

    @patch("viralunity.viralunity_consensus.validate_illumina_requirements")
    @patch("viralunity.viralunity_consensus.validate_consensus_requirements")
    @patch("viralunity.viralunity_consensus.get_samples_from_args")
    def test_validate_args_primer_scheme_not_set(
        self, mock_get_samples, mock_validate_consensus, mock_validate_illumina
    ):
        """Test validation succeeds when primer scheme is not provided."""
        mock_get_samples.return_value = {"sample1": ["file1.fastq"]}
        self.args["primer_scheme"] = None

        samples = validate_args(self.args)

        self.assertEqual(self.args["primer_scheme"], "NA")
        self.assertIn("sample1", samples)

    @patch("viralunity.viralunity_consensus.validate_illumina_requirements")
    @patch("viralunity.viralunity_consensus.validate_consensus_requirements")
    @patch("viralunity.viralunity_consensus.get_samples_from_args")
    def test_validate_args_illumina_adapters_not_exist(
        self, mock_get_samples, mock_validate_consensus, mock_validate_illumina
    ):
        """Test validation fails when adapters file doesn't exist."""
        mock_get_samples.return_value = {"sample1": ["file1.fastq"]}
        mock_validate_consensus.return_value = None
        mock_validate_illumina.side_effect = AdaptersNotFoundError(
            "Illumina adapter sequences file does not exist"
        )

        with self.assertRaises(AdaptersNotFoundError):
            validate_args(self.args)

    @patch("viralunity.viralunity_consensus.validate_illumina_requirements")
    @patch("viralunity.viralunity_consensus.validate_consensus_requirements")
    @patch("viralunity.viralunity_consensus.get_samples_from_args")
    def test_validate_args_illumina_adapters_not_set(
        self, mock_get_samples, mock_validate_consensus, mock_validate_illumina
    ):
        """Test validation fails when adapters are not provided for Illumina."""
        mock_get_samples.return_value = {"sample1": ["file1.fastq"]}
        mock_validate_consensus.return_value = None
        mock_validate_illumina.side_effect = AdaptersNotFoundError(
            "Illumina adapter sequences file is required"
        )

        self.args["adapters"] = None

        with self.assertRaises(AdaptersNotFoundError):
            validate_args(self.args)


class Test_ValidateSampleSheet(unittest.TestCase):
    """Tests for validate_sample_sheet function.

    Note: This function is now in validators module, but we test it here
    for backward compatibility with existing tests.
    """

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
        """Test sample sheet validation for Illumina data."""
        from viralunity.validators import validate_sample_sheet

        with patch("viralunity.validators.pd.read_csv", return_value=self.samplesheet_dataframe):
            samples = validate_sample_sheet("sample_sheet.csv", "illumina")
            self.assertEqual(
                samples,
                {
                    "sample1": ["R1_sample1.fastq", "R2_sample1.fastq"],
                    "sample2": ["R1_sample2.fastq", "R2_sample2.fastq"],
                },
            )

    @patch("os.path.isfile", return_value=False)
    def test_not_validate_sample_sheet_illumina(self, mock_isfile):
        """Test sample sheet validation fails when file doesn't exist."""
        from viralunity.validators import validate_sample_sheet

        with patch("viralunity.validators.pd.read_csv", return_value=self.samplesheet_dataframe):
            with self.assertRaises(Exception):
                validate_sample_sheet("sample_sheet.csv", "illumina")

    @patch("os.path.isfile", return_value=True)
    def test_validate_sample_sheet_nanopore(self, mock_isfile):
        """Test sample sheet validation for Nanopore data."""
        from viralunity.validators import validate_sample_sheet

        with patch("viralunity.validators.pd.read_csv", return_value=self.samplesheet_dataframe):
            samples = validate_sample_sheet("sample_sheet.csv", "nanopore")
            self.assertEqual(
                samples,
                {
                    "sample1": ["R1_sample1.fastq"],
                    "sample2": ["R1_sample2.fastq"],
                },
            )

    @patch("os.path.isfile", return_value=False)
    def test_not_validate_sample_sheet_nanopore(self, mock_isfile):
        """Test sample sheet validation fails when file doesn't exist."""
        from viralunity.validators import validate_sample_sheet

        with patch("viralunity.validators.pd.read_csv", return_value=self.samplesheet_dataframe):
            with self.assertRaises(Exception):
                validate_sample_sheet("sample_sheet.csv", "nanopore")


class Test_GenerateConfigFile(unittest.TestCase):
    def setUp(self):
        self.args = {
            "sample_sheet": "sample_sheet.csv",
            "config_file": "config_file.yaml",
            "output": "output_dir",
            "run_name": "run_name",
            "reference": "reference.fasta",
            "segmented_reference": None,
            "primer_scheme": "scheme",
            "minimum_coverage": 5,
            "adapters": "adapters.fasta",
            "minimum_read_length": 50,
            "trim": 0,
            "create_config_only": False,
            "threads": 1,
            "threads_total": 1,
        }

    @patch("builtins.open", new_callable=mock_open)
    @patch("os.makedirs")
    @patch("viralunity.config_generator.yaml.dump")
    def test_generate_config_file_illumina(self, mock_yaml_dump, mock_makedirs, mock_open):
        """Test config file generation for Illumina data."""
        self.args["data_type"] = "illumina"
        self.samples = {
            "sample1": ["R1_sample1.fastq", "R2_sample1.fastq"],
            "sample2": ["R1_sample2.fastq", "R2_sample2.fastq"],
        }

        generate_config_file(self.samples, self.args)

        # config_file.yaml has no directory component, so makedirs should NOT be called
        mock_makedirs.assert_not_called()
        mock_open.assert_called_once_with("config_file.yaml", "w")
        # Check that yaml.dump was called (once per section)
        self.assertGreaterEqual(mock_yaml_dump.call_count, 1)
        # Aggregate all dumped sections into one dict
        config_dict = {}
        for call in mock_yaml_dump.call_args_list:
            config_dict.update(call[0][0])
        self.assertIn("samples", config_dict)
        self.assertEqual(config_dict["data"], "illumina")
        self.assertEqual(config_dict["reference"], "reference.fasta")
        self.assertEqual(config_dict["scheme"], "scheme")
        self.assertEqual(config_dict["minimum_depth"], 5)
        self.assertEqual(config_dict["threads"], 1)
        self.assertEqual(config_dict["workflow_path"], sys.path[0])
        self.assertEqual(config_dict["output"], "output_dir/run_name/")
        self.assertEqual(config_dict["adapters"], "adapters.fasta")
        self.assertEqual(config_dict["minimum_length"], 50)
        self.assertEqual(config_dict["trim_head"], 0)
        self.assertEqual(config_dict["trim_tail"], 0)
        self.assertEqual(config_dict["run_isnv"], False)

    @patch("builtins.open", new_callable=mock_open)
    @patch("os.makedirs")
    @patch("viralunity.config_generator.yaml.dump")
    def test_generate_config_file_illumina_with_isnv(
        self, mock_yaml_dump, mock_makedirs, mock_open
    ):
        """Test config file generation for Illumina data with iSNV enabled."""
        self.args["data_type"] = "illumina"
        self.args["run_isnv"] = True
        self.samples = {
            "sample1": ["R1_sample1.fastq", "R2_sample1.fastq"],
            "sample2": ["R1_sample2.fastq", "R2_sample2.fastq"],
        }

        generate_config_file(self.samples, self.args)

        self.assertGreaterEqual(mock_yaml_dump.call_count, 1)
        config_dict = {}
        for call in mock_yaml_dump.call_args_list:
            config_dict.update(call[0][0])
        self.assertEqual(config_dict["data"], "illumina")
        self.assertEqual(config_dict["run_isnv"], True)

    @patch("builtins.open", new_callable=mock_open)
    @patch("os.makedirs")
    @patch("viralunity.config_generator.yaml.dump")
    def test_generate_config_file_nanopore(self, mock_yaml_dump, mock_makedirs, mock_open):
        """Test config file generation for Nanopore data."""
        self.args["data_type"] = "nanopore"
        self.samples = {
            "sample1": ["R1_sample1.fastq"],
            "sample2": ["R1_sample2.fastq"],
        }

        generate_config_file(self.samples, self.args)

        # config_file.yaml has no directory component, so makedirs should NOT be called
        mock_makedirs.assert_not_called()
        mock_open.assert_called_once_with("config_file.yaml", "w")
        # Check that yaml.dump was called (once per section)
        self.assertGreaterEqual(mock_yaml_dump.call_count, 1)
        # Aggregate all dumped sections into one dict
        config_dict = {}
        for call in mock_yaml_dump.call_args_list:
            config_dict.update(call[0][0])
        self.assertIn("samples", config_dict)
        self.assertEqual(config_dict["data"], "nanopore")
        self.assertEqual(config_dict["reference"], "reference.fasta")
        self.assertEqual(config_dict["scheme"], "scheme")
        self.assertEqual(config_dict["minimum_depth"], 5)
        self.assertEqual(config_dict["threads"], 1)
        self.assertEqual(config_dict["workflow_path"], sys.path[0])
        self.assertEqual(config_dict["output"], "output_dir/run_name/")
        # Nanopore should not have Illumina-specific settings
        self.assertNotIn("adapters", config_dict)
        self.assertNotIn("trim_head", config_dict)
        self.assertNotIn("trim_tail", config_dict)
        # Nanopore should have nanopore-specific settings
        self.assertIn("minimum_length", config_dict)
        self.assertIn("af_threshold", config_dict)
        self.assertIn("chunk_size", config_dict)
        self.assertIn("clair3_model", config_dict)

    @patch("builtins.open", new_callable=mock_open)
    @patch("os.makedirs")
    @patch("viralunity.config_generator.yaml.dump")
    def test_generate_config_file_segmented_reference(
        self, mock_yaml_dump, mock_makedirs, mock_open
    ):
        """Test config file generation with segmented reference."""
        self.args["data_type"] = "nanopore"
        self.args["reference"] = {"S": "/path/to/S.fasta", "L": "/path/to/L.fasta"}
        self.samples = {
            "sample1": ["R1_sample1.fastq"],
        }

        generate_config_file(self.samples, self.args)

        self.assertGreaterEqual(mock_yaml_dump.call_count, 1)
        config_dict = {}
        for call in mock_yaml_dump.call_args_list:
            config_dict.update(call[0][0])
        self.assertEqual(
            config_dict["reference"], {"S": "/path/to/S.fasta", "L": "/path/to/L.fasta"}
        )


class Test_MainFunction(unittest.TestCase):
    @patch(
        "viralunity.viralunity_consensus.validate_args",
        return_value={"sample1": ["file1.fastq"]},
    )
    @patch("viralunity.viralunity_consensus.generate_config_file")
    @patch("viralunity.viralunity_consensus.run_snakemake_workflow", return_value=True)
    def test_main_success(self, mock_run_workflow, mock_generate_config_file, mock_validate_args):
        """Test main function succeeds when workflow completes."""
        result = main(
            {
                "config_file": "config_file.yaml",
                "threads_total": 1,
                "data_type": "illumina",
                "create_config_only": False,
            }
        )
        self.assertEqual(result, 0)
        mock_run_workflow.assert_called_once()

    @patch(
        "viralunity.viralunity_consensus.validate_args",
        return_value={"sample1": ["file1.fastq"]},
    )
    @patch("viralunity.viralunity_consensus.generate_config_file")
    @patch("viralunity.viralunity_consensus.run_snakemake_workflow", return_value=False)
    def test_main_create_config_only(
        self, mock_run_workflow, mock_generate_config_file, mock_validate_args
    ):
        """Test main function exits early when create_config_only is True."""
        result = main(
            {
                "config_file": "config_file.yaml",
                "threads_total": 1,
                "data_type": "illumina",
                "create_config_only": True,
            }
        )
        self.assertEqual(result, 0)
        mock_run_workflow.assert_not_called()

    @patch(
        "viralunity.viralunity_consensus.validate_args",
        return_value={"sample1": ["file1.fastq"]},
    )
    @patch("viralunity.viralunity_consensus.generate_config_file")
    @patch("viralunity.viralunity_consensus.run_snakemake_workflow", return_value=False)
    def test_main_workflow_failure(
        self, mock_run_workflow, mock_generate_config_file, mock_validate_args
    ):
        """Test main function returns error code when workflow fails."""
        result = main(
            {
                "config_file": "config_file.yaml",
                "threads_total": 1,
                "data_type": "illumina",
                "create_config_only": False,
            }
        )
        self.assertEqual(result, 1)
        mock_run_workflow.assert_called_once()


if __name__ == "__main__":
    unittest.main()
