import unittest
from unittest.mock import mock_open, patch
import pandas
from io import StringIO
from viralunity.mutation_mapper import (
    get_args,
    validate_args,
    get_annotation_info,
    read_and_annotate_sequences,
    compare_sequences,
)


class Test_GetArgs(unittest.TestCase):

    def test_required_args_fail_when_only_optional_set(self):
        args = [
            "--flag",
            "CDS",
        ]
        with self.assertRaises(SystemExit):
            get_args(args)

    def test_required_args_success_when_only_required_set(self):
        args = [
            "--input",
            "input.fasta",
            "--output",
            "output.txt",
            "--annotation-file",
            "annotation.gff",
        ]
        readArgs = get_args(args)
        self.assertDictContainsSubset(
            {
                "input": "input.fasta",
                "output": "output.txt",
                "annotation_file": "annotation.gff",
            },
            readArgs,
        )

    def test_default_values_optional_args(self):
        args = [
            "--input",
            "input.fasta",
            "--output",
            "output.txt",
            "--annotation-file",
            "annotation.gff",
        ]
        readArgs = get_args(args)
        self.assertDictEqual(
            readArgs,
            {
                "input": "input.fasta",
                "output": "output.txt",
                "annotation_file": "annotation.gff",
                "flag": "gene",
            },
            readArgs,
        )


class Test_ValidateArgs(unittest.TestCase):
    def setUp(self):
        self.args = {
            "input": "input.fasta",
            "output": "output.txt",
            "annotation_file": "annotation.gff",
            "flag": "gene",
        }

    @patch("os.path.isfile", side_effect=[True, False, True])
    def test_validate_args_success(self, mock_isfile):
        samples = validate_args(self.args)
        self.assertEqual(samples, None)

    @patch("os.path.isfile", side_effect=[False, False, True])
    def test_validate_args_input_not_exist(self, mock_isfile):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, True, True])
    def test_validate_args_outout_already_exist(self, mock_isfile):
        with self.assertRaises(SystemExit):
            validate_args(self.args)

    @patch("os.path.isfile", side_effect=[True, False, False])
    def test_validate_args_annotation_file_not_exist(self, mock_isfile):
        with self.assertRaises(SystemExit):
            validate_args(self.args)


class Test_GetAnnotationInfo(unittest.TestCase):
    def setUp(self):
        self.annotation_content = """##gff-version 3
chr1\t.\tgene\t1300\t1500\t.\t+\t.\tID=gene1;Name=Gene1;
chr1\t.\ttranscript\t1050\t1500\t.\t+\t.\tID=gene2;Name=Gene2;
chr1\t.\tgene\t3000\t3902\t.\t+\t.\tID=gene3;Name=Gene3;
"""

        self.args = {
            "input": "input.fasta",
            "output": "output.txt",
            "annotation_file": StringIO(self.annotation_content),
            "flag": "gene",
        }

    def test_get_annotation_info(self):

        annotation_info = get_annotation_info(self.args)
        expected_info = {
            "Gene1": (1300, 1500),
            "Gene3": (3000, 3902),
        }
        self.assertDictEqual(annotation_info, expected_info)

# TODO: Understand with Filipe how to correctly test this function
class Test_ReadAndAnnotateSequences(unittest.TestCase):
    def setUp(self):
        self.sequences = [
            {"id": "ref1", "seq": "ATCGCGAGGCTTAACGATCG"},
            {"id": "ref2", "seq": "ATCGNCAGGCTTXXCG-TCA"},
        ]

    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data=""">ref_seq1
ATCGCGAGGCTTAACGATCG
>ref_seq2
ATCGNCAGGCTTXXCG-TCA
>ref_seq3
ATTCGGCGATAGCGGATN-X
""",
    )
    def test_read_and_annotate_sequences(self, mock_open):

        genomic_mismatches = read_and_annotate_sequences(
            {"input": "seq.fasta"},
            {"ref_seq1": (1, 20), "ref_seq2": (1, 15), "ref_seq3": (5, 20)},
        )
        self.assertEqual(genomic_mismatches, {"seq": ["G5C", "G19A"]})

class Test_CompareSequences(unittest.TestCase):
    def test_compare_sequences(self):
        ref_seq = "ATCGCGAGGCTTAACGATCG"
        seq = "ATCGNCAGGCTTXXCG-TCA"
        mutations = compare_sequences(ref_seq, seq)
        self.assertEqual(mutations, ["G5C", "G19A"])


if __name__ == "__main__":
    unittest.main()
