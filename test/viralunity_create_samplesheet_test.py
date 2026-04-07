"""Tests for viralunity create-samplesheet CLI command (click-based)."""

import unittest
from unittest.mock import patch
from click.testing import CliRunner
from viralunity.viralunity_create_samplesheet import (
    create_samplesheet,
    validate_args,
    generate_sample_sheet,
)


class Test_CreateSamplesheeetCommand(unittest.TestCase):
    """Tests for `viralunity create-samplesheet`."""

    def setUp(self):
        self.runner = CliRunner()

    def _invoke(self, extra_args=None):
        args = ["--input", "input/dir", "--output", "output.file"] + (extra_args or [])
        with patch("viralunity.viralunity_create_samplesheet.validate_args"), patch(
            "viralunity.viralunity_create_samplesheet.generate_sample_sheet"
        ):
            return self.runner.invoke(create_samplesheet, args, catch_exceptions=False)

    def test_get_args_required(self):
        """Missing required args (--input, --output) should exit non-zero."""
        result = self.runner.invoke(create_samplesheet, ["--level", "1"])
        self.assertNotEqual(result.exit_code, 0)

    def test_required_args_success_when_only_required_set(self):
        result = self._invoke()
        self.assertEqual(result.exit_code, 0, result.output)

    def test_default_values_optional_args(self):
        """Check that defaults match the expected values."""
        with patch(
            "viralunity.viralunity_create_samplesheet.validate_args"
        ) as mock_validate, patch(
            "viralunity.viralunity_create_samplesheet.generate_sample_sheet"
        ) as mock_gen:
            result = self.runner.invoke(
                create_samplesheet,
                ["--input", "input/dir", "--output", "output.file"],
                catch_exceptions=False,
            )
        self.assertEqual(result.exit_code, 0, result.output)
        called_args = mock_validate.call_args[0][0]
        self.assertEqual(called_args["input"], "input/dir")
        self.assertEqual(called_args["output"], "output.file")
        self.assertEqual(called_args["level"], 1)
        self.assertEqual(called_args["pattern"], "R1")
        self.assertEqual(called_args["separator"], "-")


class Test_ValidateArgs(unittest.TestCase):
    def setUp(self):
        self.args = {
            "input": "input/dir",
            "output": "output.file",
            "level": 1,
        }

    @patch("os.path.isdir", return_value=True)
    @patch("os.path.isfile", return_value=False)
    def test_validate_args_success(self, mock_isfile, mock_isdir):
        validated_args = validate_args(self.args)
        self.assertEqual(validated_args, None)

    @patch("os.path.isdir", return_value=False)
    def test_validate_args_input_not_exist(self, mock_isdir):
        with self.assertRaises(Exception):
            validate_args(self.args)

    @patch("os.path.isdir", return_value=True)
    @patch("os.path.isfile", return_value=True)
    def test_validate_args_output_exists(self, mock_isfile, mock_isdir):
        with self.assertRaises(Exception):
            validate_args(self.args)


class Test_GenerateSamplesheet(unittest.TestCase):
    def setUp(self):
        self.args = {
            "input": "input/dir",
            "output": "output.file",
            "separator": "_",
            "pattern": "R1",
        }

    @patch(
        "glob.glob",
        side_effect=[
            ["input/dir/1", "input/dir/2"],
            ["input/dir/1/R1_sample1.fastq", "input/dir/1/R1_sample2.fastq"],
            ["input/dir/2/R2_sample1.fastq", "input/dir/2/R2_sample2.fastq"],
        ],
    )
    @patch("os.path.isdir", return_value=True)
    @patch("os.path.isfile", return_value=True)
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    def test_generate_samplesheet_level_1(
        self, mock_open, mock_isfile, mock_isdir, mock_glob
    ):
        self.args["level"] = 1
        generate_sample_sheet(self.args)
        mock_open.assert_called_with("output.file", "w")
        handle = mock_open()
        handle.write.assert_any_call(
            "1,input/dir/1/R1_sample1.fastq,input/dir/1/R1_sample2.fastq\n"
        )
        handle.write.assert_any_call(
            "2,input/dir/2/R2_sample1.fastq,input/dir/2/R2_sample2.fastq\n"
        )

    @patch(
        "glob.glob",
        side_effect=[
            ["input/dir/1/R1_sample1.fastq", "input/dir/1/R1_sample2.fastq"],
            ["input/dir/1/R1_sample1.fastq", "input/dir/1/R1_sample2.fastq"],
            ["input/dir/1/R1_sample1.fastq", "input/dir/1/R1_sample2.fastq"],
        ],
    )
    @patch("os.path.isfile", return_value=True)
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    def test_generate_samplesheet_level_0(self, mock_open, mock_isfile, mock_glob):
        self.args["level"] = 0
        generate_sample_sheet(self.args)
        mock_open.assert_called_with("output.file", "w")
        handle = mock_open()
        handle.write.assert_any_call(
            "R1,input/dir/1/R1_sample1.fastq,input/dir/1/R1_sample2.fastq\n"
        )


if __name__ == "__main__":
    unittest.main()
