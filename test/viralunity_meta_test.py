import unittest
from unittest.mock import patch, mock_open
from viralunity.viralunity_meta import (
    validate_args,
    generate_config_file,
    main,
)
from viralunity.exceptions import (
    Kraken2DatabaseNotFoundError,
    KronaDatabaseNotFoundError,
    SampleConfigurationNotFoundError,
    AdaptersNotFoundError
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

    @patch("viralunity.viralunity_meta.validate_illumina_requirements")
    @patch("viralunity.viralunity_meta.validate_metagenomics_requirements")
    @patch("viralunity.viralunity_meta.get_samples_from_args")
    def test_validate_args_success(
        self, mock_get_samples, mock_validate_meta, mock_validate_illumina
    ):
        """Test successful validation with all validators passing."""
        mock_get_samples.return_value = {"sample1": ["file1_R1.fastq", "file1_R2.fastq"]}
        
        samples = validate_args(self.args)
        
        self.assertIn("sample1", samples)
        mock_get_samples.assert_called_once_with(self.args)
        mock_validate_meta.assert_called_once_with(self.args)
        mock_validate_illumina.assert_called_once_with(self.args)

    @patch("viralunity.viralunity_meta.get_samples_from_args")
    def test_validate_args_sample_sheet_not_exist(self, mock_get_samples):
        """Test validation fails when sample sheet cannot be retrieved."""
        mock_get_samples.side_effect = SampleConfigurationNotFoundError("Sample sheet not found")
        
        with self.assertRaises(SampleConfigurationNotFoundError):
            validate_args(self.args)
        
        mock_get_samples.assert_called_once_with(self.args)

    @patch("viralunity.viralunity_meta.validate_metagenomics_requirements")
    @patch("viralunity.viralunity_meta.get_samples_from_args")
    def test_validate_args_kraken2_db_not_exist(
        self, mock_get_samples, mock_validate_meta
    ):
        """Test validation fails when Kraken2 database doesn't exist."""
        mock_get_samples.return_value = {"sample1": ["file1.fastq"]}
        mock_validate_meta.side_effect = Kraken2DatabaseNotFoundError("Kraken2 database directory does not exist")
        
        with self.assertRaises(Kraken2DatabaseNotFoundError):
            validate_args(self.args)

    @patch("viralunity.viralunity_meta.validate_metagenomics_requirements")
    @patch("viralunity.viralunity_meta.get_samples_from_args")
    def test_validate_args_krona_db_not_exist(
        self, mock_get_samples, mock_validate_meta
    ):
        """Test validation fails when Krona database doesn't exist."""
        mock_get_samples.return_value = {"sample1": ["file1.fastq"]}
        mock_validate_meta.side_effect = KronaDatabaseNotFoundError("Krona database directory does not exist")
        
        with self.assertRaises(KronaDatabaseNotFoundError):
            validate_args(self.args)

    @patch("viralunity.viralunity_meta.validate_illumina_requirements")
    @patch("viralunity.viralunity_meta.validate_metagenomics_requirements")
    @patch("viralunity.viralunity_meta.get_samples_from_args")
    def test_validate_args_adapters_not_exist(
        self, mock_get_samples, mock_validate_meta, mock_validate_illumina
    ):
        """Test validation fails when adapters file doesn't exist."""
        mock_get_samples.return_value = {"sample1": ["file1.fastq"]}
        mock_validate_meta.return_value = None
        mock_validate_illumina.side_effect = AdaptersNotFoundError("Illumina adapter sequences file does not exist")
        
        with self.assertRaises(AdaptersNotFoundError):
            validate_args(self.args)

    @patch("viralunity.viralunity_meta.validate_illumina_requirements")
    @patch("viralunity.viralunity_meta.validate_metagenomics_requirements")
    @patch("viralunity.viralunity_meta.get_samples_from_args")
    def test_validate_args_config_file_exists(
        self, mock_get_samples, mock_validate_meta, mock_validate_illumina
    ):
        """Test validation succeeds even if config file already exists."""
        mock_get_samples.return_value = {"sample1": ["file1_R1.fastq", "file1_R2.fastq"]}
        
        samples = validate_args(self.args)
        
        self.assertIn("sample1", samples)

    @patch("viralunity.viralunity_meta.validate_illumina_requirements")
    @patch("viralunity.viralunity_meta.validate_metagenomics_requirements")
    @patch("viralunity.viralunity_meta.get_samples_from_args")
    def test_validate_args_output_dir_exists(
        self, mock_get_samples, mock_validate_meta, mock_validate_illumina
    ):
        """Test validation succeeds even if output directory already exists."""
        mock_get_samples.return_value = {"sample1": ["file1_R1.fastq", "file1_R2.fastq"]}
        
        samples = validate_args(self.args)
        
        self.assertIn("sample1", samples)


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

    @patch("builtins.open", new_callable=mock_open)
    @patch("os.makedirs")
    @patch("viralunity.config_generator.yaml.dump")
    def test_generate_config_file_illumina(self, mock_yaml_dump, mock_makedirs, mock_open):
        self.args["data_type"] = "illumina"
        self.samples = {
            "sample1": ["sample1_R1.fastq", "sample1_R2.fastq"],
            "sample2": ["sample2_R1.fastq", "sample2_R2.fastq"],
        }
        generate_config_file(self.samples, self.args)
        # makedirs is only called if config file path has a directory component
        config_dir = os.path.dirname("config_file.yaml")
        if config_dir:
            mock_makedirs.assert_called_once_with(config_dir, exist_ok=True)
        mock_open.assert_called_once_with("config_file.yaml", "w")
        # Check that yaml.dump was called with correct config structure
        self.assertEqual(mock_yaml_dump.call_count, 1)
        call_args = mock_yaml_dump.call_args
        config_dict = call_args[0][0]  # First positional argument
        self.assertIn("samples", config_dict)
        self.assertEqual(config_dict["data"], "illumina")
        self.assertEqual(config_dict["kraken2_database"], "kraken2_db")
        self.assertEqual(config_dict["krona_database"], "krona_db")
        self.assertEqual(config_dict["threads"], 1)
        self.assertEqual(config_dict["output"], "output_dir/run_name/")
        self.assertEqual(config_dict["adapters"], "adapters.fasta")
        self.assertEqual(config_dict["minimum_length"], 50)
        self.assertEqual(config_dict["trim"], 0)
        self.assertEqual(config_dict["remove_human_reads"], True)
        self.assertEqual(config_dict["remove_unclassified_reads"], False)

    @patch("builtins.open", new_callable=mock_open)
    @patch("os.makedirs")
    @patch("viralunity.config_generator.yaml.dump")
    def test_generate_config_file_nanopore(self, mock_yaml_dump, mock_makedirs, mock_open):
        self.args["data_type"] = "nanopore"
        self.samples = {
            "sample1": ["sample1.fastq"],
            "sample2": ["sample2.fastq"],
        }
        generate_config_file(self.samples, self.args)
        # makedirs is only called if config file path has a directory component
        config_dir = os.path.dirname("config_file.yaml")
        if config_dir:
            mock_makedirs.assert_called_once_with(config_dir, exist_ok=True)
        mock_open.assert_called_once_with("config_file.yaml", "w")
        # Check that yaml.dump was called with correct config structure
        self.assertEqual(mock_yaml_dump.call_count, 1)
        call_args = mock_yaml_dump.call_args
        config_dict = call_args[0][0]  # First positional argument
        self.assertIn("samples", config_dict)
        self.assertEqual(config_dict["data"], "nanopore")
        self.assertEqual(config_dict["kraken2_database"], "kraken2_db")
        self.assertEqual(config_dict["krona_database"], "krona_db")
        self.assertEqual(config_dict["threads"], 1)
        self.assertEqual(config_dict["output"], "output_dir/run_name/")
        self.assertEqual(config_dict["remove_human_reads"], True)
        self.assertEqual(config_dict["remove_unclassified_reads"], False)
        # Nanopore should not have Illumina-specific settings
        self.assertNotIn("adapters", config_dict)
        self.assertNotIn("minimum_length", config_dict)
        self.assertNotIn("trim", config_dict)


class Test_MainFunction(unittest.TestCase):
    @patch("viralunity.viralunity_meta.validate_args", return_value={"sample1": ["file1.fastq"]})
    @patch("viralunity.viralunity_meta.generate_config_file")
    @patch("viralunity.viralunity_meta.run_snakemake_workflow", return_value=True)
    def test_main_success(
        self,
        mock_run_workflow,
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
        mock_run_workflow.assert_called_once()

    @patch("viralunity.viralunity_meta.validate_args", return_value={"sample1": ["file1.fastq"]})
    @patch("viralunity.viralunity_meta.generate_config_file")
    @patch("viralunity.viralunity_meta.run_snakemake_workflow", return_value=True)
    def test_main_create_config_only(
        self,
        mock_run_workflow,
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
        mock_run_workflow.assert_not_called()

    @patch("viralunity.viralunity_meta.validate_args", return_value={"sample1": ["file1.fastq"]})
    @patch("viralunity.viralunity_meta.generate_config_file")
    @patch("viralunity.viralunity_meta.run_snakemake_workflow", return_value=False)
    def test_main_failure(
        self,
        mock_run_workflow,
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
        mock_run_workflow.assert_called_once()


if __name__ == "__main__":
    unittest.main()
