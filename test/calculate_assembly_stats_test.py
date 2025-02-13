import unittest
import subprocess
import pandas as pd
from unittest.mock import patch
from viralunity.calculate_assembly_stats import (
    get_number_of_reads,
    get_number_of_mapped_reads,
    get_coverage_info,
    generate_output,
    main,
)


class TestGetNumberOfReads(unittest.TestCase):

    @patch("subprocess.Popen")
    def test_get_number_of_reads_gz(self, mock_popen):
        mock_popen.return_value.communicate.return_value = (b"4\n", None)
        fastq = "sample.fastq.gz"
        result = get_number_of_reads(fastq)
        mock_popen.assert_called_with(
            'gunzip -c sample.fastq.gz | grep -cE "^\+$"',
            stdout=subprocess.PIPE,
            shell=True,
        )
        self.assertEqual(result, 4)

    @patch("subprocess.Popen")
    def test_get_number_of_reads(self, mock_popen):
        mock_popen.return_value.communicate.return_value = (b"4\n", None)
        fastq = "sample.fastq"
        result = get_number_of_reads(fastq)
        mock_popen.assert_called_with(
            'grep -cE "^\+$" sample.fastq', stdout=subprocess.PIPE, shell=True
        )
        self.assertEqual(result, 4)


class TestGetNumberOfMappedReads(unittest.TestCase):

    @patch("subprocess.Popen")
    def test_get_number_of_mapped_reads(self, mock_popen):
        mock_popen.return_value.communicate.return_value = (b"10\n", None)
        bam = "sample.bam"
        result = get_number_of_mapped_reads(bam)
        mock_popen.assert_called_with(
            "samtools view -c -F 260 sample.bam", stdout=subprocess.PIPE, shell=True
        )
        self.assertEqual(result, 10)


class TestGetCoverageInfo(unittest.TestCase):

    @patch("pandas.read_csv")
    def test_get_coverage_info(self, mock_read_csv):
        data = {
            0: [1, 2, 3, 4, 5],
            1: [10, 20, 30, 40, 50],
            2: [0, 10, 100, 1000, 10000],
        }
        df = pd.DataFrame(data)
        mock_read_csv.return_value = df
        table_cov = "coverage.txt"
        minimum_depth = 20

        result = get_coverage_info(table_cov, minimum_depth)
        expected_average_depth = 2222
        expected_percentage_of_sites_above_10x = 0.8
        expected_percentage_of_sites_above_100x = 0.6
        expected_percentage_of_sites_above_1000x = 0.4
        expected_percentage_of_sites_above_specified_threshold = 0.6

        self.assertEqual(result[0], expected_average_depth)
        self.assertEqual(result[1], expected_percentage_of_sites_above_10x)
        self.assertEqual(result[2], expected_percentage_of_sites_above_100x)
        self.assertEqual(result[3], expected_percentage_of_sites_above_1000x)
        self.assertEqual(
            result[4], expected_percentage_of_sites_above_specified_threshold
        )


class TestWriteOutput(unittest.TestCase):

    @patch("viralunity.calculate_assembly_stats.get_number_of_reads")
    @patch("viralunity.calculate_assembly_stats.get_number_of_mapped_reads")
    @patch("viralunity.calculate_assembly_stats.get_coverage_info")
    @patch("pandas.DataFrame.to_csv")
    def test_write_output(
        self,
        mock_to_csv,
        mock_get_coverage_info,
        mock_get_number_of_mapped_reads,
        mock_get_number_of_reads,
    ):
        mock_get_number_of_reads.side_effect = [100, 80]
        mock_get_number_of_mapped_reads.return_value = 70
        mock_get_coverage_info.return_value = (30, 0.9, 0.8, 0.7, 0.6)

        fastq = "sample.fastq"
        trim_fastq = "sample_trim.fastq"
        bam = "sample.bam"
        table_cov = "coverage.txt"
        sample_name = "sample"
        minimum_depth = 20
        output = "output.csv"

        df_out = generate_output(
            fastq, trim_fastq, bam, table_cov, sample_name, minimum_depth, output
        )

        mock_get_number_of_reads.assert_any_call(fastq)
        mock_get_number_of_reads.assert_any_call(trim_fastq)
        mock_get_number_of_mapped_reads.assert_called_once_with(bam)
        mock_get_coverage_info.assert_called_once_with(table_cov, minimum_depth)

        expected_data = [
            sample_name,
            100,
            80,
            70,
            30,
            0.9,
            0.8,
            0.7,
            0.6,
        ]
        expected_df = pd.DataFrame(expected_data).T
        assert df_out.equals(expected_df)


class TestMainFunction(unittest.TestCase):

    @patch("viralunity.calculate_assembly_stats.generate_output")
    @patch("pandas.DataFrame.to_csv")
    def test_main(self, mock_to_csv, mock_generate_output):
        inputs = ["sample.fastq", "", "sample_trim.fastq", "my/dir/sample.sorted.bam", "coverage.txt"]
        output = "output.csv"
        minimum_depth = 20

        mock_generate_output.return_value = pd.DataFrame(
            [["sample", 100, 80, 70, 30, 0.9, 0.8, 0.7, 0.6]]
        )

        main(inputs, output, minimum_depth)

        mock_generate_output.assert_called_once_with(
            "sample.fastq",
            "sample_trim.fastq",
            "my/dir/sample.sorted.bam",
            "coverage.txt",
            "sample",
            minimum_depth,
            output,
        )
        mock_to_csv.assert_called_once_with(output, header=False, index=False)


if __name__ == "__main__":
    unittest.main()
