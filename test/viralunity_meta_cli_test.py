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
                "assembly_summary": None,
                "config_file": "config_file.yaml",
                "create_config_only": False,
                "data_type": "illumina",
                "diamond_database": None,
                "diamond_sensitivity": "sensitive",
                "evalue": 1e-10,
                "host_reference": "NA",
                "kraken2_database": "kraken2_db",
                "krona_database": "krona_db",
                "minimum_read_length": 50,
                "negative_controls": "",
                "negative_p_threshold": 0.01,
                "output": "output_dir",
                "pipeline": "v1",
                "remove_human_reads": False,
                "remove_human_sequences": False,
                "remove_unclassified_reads": False,
                "remove_unclassified_sequences": False,
                "run_diamond": False,
                "run_diamond_contigs": False,
                "run_diamond_reads": False,
                "run_denovo_assembly": False,
                "run_kraken2": False,
                "run_kraken2_contigs": False,
                "run_kraken2_reads": False,
                "run_name": "undefined",
                "sample_sheet": "sample_sheet.csv",
                "taxdump": None,
                "taxid_to_family": None,
                "threads": 1,
                "threads_total": 1,
                "trim": 0,
                "bleed_fraction": 0.005,
            },
        )

    def test_v2_args_default_values(self):
        """Test that v2 arguments have correct default values when not specified."""
        args = [
            "meta",
            "--data-type",
            "nanopore",
            "--pipeline",
            "v2",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--krona-database",
            "krona_db",
        ]
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_meta(subparsers)
        readArgs = vars(parser.parse_args(args))
        
        # Check v2 workflow toggles default to False
        self.assertFalse(readArgs.get("run_denovo_assembly", True))
        self.assertFalse(readArgs.get("run_kraken2_contigs", True))
        self.assertFalse(readArgs.get("run_diamond_contigs", True))
        self.assertFalse(readArgs.get("run_kraken2_reads", True))
        self.assertFalse(readArgs.get("run_diamond_reads", True))
        self.assertFalse(readArgs.get("run_kraken2", True))
        self.assertFalse(readArgs.get("run_diamond", True))
        
        # Check DIAMOND parameters defaults
        self.assertIsNone(readArgs.get("diamond_database"))
        self.assertEqual(readArgs.get("diamond_sensitivity"), "sensitive")
        self.assertEqual(readArgs.get("evalue"), 1e-10)
        
        # Check host removal default
        self.assertEqual(readArgs.get("host_reference"), "NA")
        
        # Check taxonomy resources default to None
        self.assertIsNone(readArgs.get("taxdump"))
        self.assertIsNone(readArgs.get("assembly_summary"))
        self.assertIsNone(readArgs.get("taxid_to_family"))

    def test_v2_args_when_explicitly_set(self):
        """Test that v2 arguments are correctly parsed when explicitly provided."""
        args = [
            "meta",
            "--data-type",
            "nanopore",
            "--pipeline",
            "v2",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--krona-database",
            "krona_db",
            "--run-denovo-assembly",
            "--run-kraken2-contigs",
            "--run-diamond-contigs",
            "--run-kraken2-reads",
            "--run-diamond-reads",
            "--diamond-database",
            "/path/to/diamond.db",
            "--diamond-sensitivity",
            "very-sensitive",
            "--evalue",
            "1e-5",
            "--host-reference",
            "/path/to/host.fasta",
            "--taxdump",
            "/path/to/taxdump",
            "--assembly-summary",
            "/path/to/assembly_summary.tsv",
            "--taxid-to-family",
            "/path/to/taxid_to_family.csv",
        ]
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_meta(subparsers)
        readArgs = vars(parser.parse_args(args))
        
        # Check v2 workflow toggles are True when set
        self.assertTrue(readArgs.get("run_denovo_assembly"))
        self.assertTrue(readArgs.get("run_kraken2_contigs"))
        self.assertTrue(readArgs.get("run_diamond_contigs"))
        self.assertTrue(readArgs.get("run_kraken2_reads"))
        self.assertTrue(readArgs.get("run_diamond_reads"))
        
        # Check DIAMOND parameters
        self.assertEqual(readArgs.get("diamond_database"), "/path/to/diamond.db")
        self.assertEqual(readArgs.get("diamond_sensitivity"), "very-sensitive")
        self.assertEqual(readArgs.get("evalue"), 1e-5)
        
        # Check host removal
        self.assertEqual(readArgs.get("host_reference"), "/path/to/host.fasta")
        
        # Check taxonomy resources
        self.assertEqual(readArgs.get("taxdump"), "/path/to/taxdump")
        self.assertEqual(readArgs.get("assembly_summary"), "/path/to/assembly_summary.tsv")
        self.assertEqual(readArgs.get("taxid_to_family"), "/path/to/taxid_to_family.csv")

    def test_v2_args_present_with_v1_pipeline(self):
        """Test that v2 arguments are available even when using v1 pipeline (they're optional)."""
        args = [
            "meta",
            "--data-type",
            "illumina",
            "--pipeline",
            "v1",
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
            "--run-diamond-reads",
            "--diamond-database",
            "/path/to/diamond.db",
        ]
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_meta(subparsers)
        readArgs = vars(parser.parse_args(args))
        
        # Verify v2 args are accessible even with v1 pipeline
        self.assertTrue(readArgs.get("run_diamond_reads"))
        self.assertEqual(readArgs.get("diamond_database"), "/path/to/diamond.db")
        self.assertEqual(readArgs.get("pipeline"), "v1")

    def test_v2_hidden_backward_compatibility_flags(self):
        """Test that hidden backward-compatibility flags work correctly."""
        args = [
            "meta",
            "--data-type",
            "nanopore",
            "--pipeline",
            "v2",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--krona-database",
            "krona_db",
            "--run-kraken2",
            "--run-diamond",
        ]
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_meta(subparsers)
        readArgs = vars(parser.parse_args(args))
        
        # Check hidden flags are set
        self.assertTrue(readArgs.get("run_kraken2"))
        self.assertTrue(readArgs.get("run_diamond"))

    def test_v2_args_with_all_common_args(self):
        """Test v2 arguments work correctly with all common arguments set."""
        args = [
            "meta",
            "--data-type",
            "nanopore",
            "--pipeline",
            "v2",
            "--sample-sheet",
            "sample_sheet.csv",
            "--config-file",
            "config_file.yaml",
            "--output",
            "output_dir",
            "--krona-database",
            "krona_db",
            "--run-name",
            "test_run",
            "--threads",
            "4",
            "--threads-total",
            "8",
            "--run-denovo-assembly",
            "--diamond-sensitivity",
            "sensitive",
            "--evalue",
            "1e-8",
            "--host-reference",
            "/path/to/host.fasta",
        ]
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        fill_arg_parser_meta(subparsers)
        readArgs = vars(parser.parse_args(args))
        
        # Verify common args
        self.assertEqual(readArgs.get("data_type"), "nanopore")
        self.assertEqual(readArgs.get("pipeline"), "v2")
        self.assertEqual(readArgs.get("run_name"), "test_run")
        self.assertEqual(readArgs.get("threads"), 4)
        self.assertEqual(readArgs.get("threads_total"), 8)
        
        # Verify v2 args
        self.assertTrue(readArgs.get("run_denovo_assembly"))
        self.assertEqual(readArgs.get("diamond_sensitivity"), "sensitive")
        self.assertEqual(readArgs.get("evalue"), 1e-8)
        self.assertEqual(readArgs.get("host_reference"), "/path/to/host.fasta")

if __name__ == "__main__":
    unittest.main()
