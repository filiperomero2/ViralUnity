import unittest
from unittest.mock import mock_open, patch
from viralunity.rename_sequences import rename_sequences

class TestRenameSequences(unittest.TestCase):
    @patch("builtins.open", new_callable=mock_open, read_data=">header\nATCG\nCGTA\nCCGC\nAACG\n")
    def test_rename_sequences(self, mock_open):
        input = "assembly/consensus/final_consensus/1234.consensus.fasta"
        output = "assembly/consensus/final_consensus/1234.consensus.renamed.fasta"
        
        rename_sequences(input, output)
        
        mock_open.assert_any_call("assembly/consensus/final_consensus/1234.consensus.fasta")
        mock_open.assert_any_call("assembly/consensus/final_consensus/1234.consensus.renamed.fasta", "w")
        handle = mock_open()
        handle.write.assert_called_once_with(">1234\nATCG\nCGTA\nCCGC\nAACG\n")

if __name__ == "__main__": # Check if this brake snakemake
    unittest.main()