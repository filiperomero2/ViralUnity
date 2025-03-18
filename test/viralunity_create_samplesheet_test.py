import unittest
from unittest.mock import patch
from viralunity.viralunity_create_samplesheet import (
    get_args,
    validate_args,
    generate_sample_sheet,
    main
)


class Test_GetArgs(unittest.TestCase):

    def test_get_args_required(self):
        args = [
            "--level",
            "1",
        ]
        with self.assertRaises(SystemExit):
            get_args(args)

    def test_required_args_success_when_only_required_set(self):
        args = ["--input", "input/dir", "--output", "output.file"]
        readArgs = get_args(args)
        self.assertDictContainsSubset(
            {
                "input": "input/dir",
                "output": "output.file",
            },
            readArgs,
        )

    def test_default_values_optional_args(self):
        args = ["--input", "input/dir", "--output", "output.file"]
        readArgs = get_args(args)
        self.assertDictEqual(
            {
                "input": "input/dir",
                "output": "output.file",
                "level": 1,
            },
            readArgs,
        )


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
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isdir", return_value=True)
    @patch("os.path.isfile", return_value=True)
    def test_validate_args_output_exists(self, mock_isfile, mock_isdir):
        with self.assertRaises(SystemExit):
            validate_args(self.args)


class Test_GenerateSamplesheet(unittest.TestCase):
    def setUp(self):
        self.args = {
            "input": "input/dir",
            "output": "output.file"
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
    def test_generate_samplesheet_level_1(self,mock_open,mock_isfile, mock_isdir, mock_glob):
        self.args['level'] = 1
        generate_sample_sheet(self.args)
        mock_open.assert_called_with("output.file", "w")
        handle = mock_open()
        handle.write.assert_any_call("1,input/dir/1/R1_sample1.fastq,input/dir/1/R1_sample2.fastq\n")
        handle.write.assert_any_call("2,input/dir/2/R2_sample1.fastq,input/dir/2/R2_sample2.fastq\n")

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
    def test_generate_samplesheet_level_0(self,mock_open,mock_isfile,mock_glob):
        self.args['level'] = 0
        generate_sample_sheet(self.args)
        mock_open.assert_called_with("output.file", "w")
        handle = mock_open()
        handle.write.assert_any_call("R1,input/dir/1/R1_sample1.fastq,input/dir/1/R1_sample2.fastq\n")

@patch("sys.argv", retuen_value=["script_name", "--input", "input/dir", "--output", "output.file"])
@patch("viralunity.create_viralunity_samplesheet.get_args", return_value={"input": "input/dir", "output": "output.file"})
@patch("viralunity.create_viralunity_samplesheet.validate_args")
@patch("viralunity.create_viralunity_samplesheet.generate_sample_sheet")
class Test_MainFunction(unittest.TestCase):
    
    def test_main_function(self, mock_generate, mock_validate, mock_get_args, mock_sys):
        main()
        mock_get_args.assert_called_once()
        mock_validate.assert_called_once()
        mock_generate.assert_called_once()
        

if __name__ == "__main__":
    unittest.main()
