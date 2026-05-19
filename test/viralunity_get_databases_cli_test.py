"""Tests for the viralunity get-databases CLI command group.

Focus on the DNA-contamination filter that protects ``diamond makedb`` from
the NCBI Datasets bug where some virus protein FASTA records are actually
nucleotide CDS sequences.
"""

import os
import tempfile
import textwrap
import unittest
from pathlib import Path

from click.testing import CliRunner

from viralunity.viralunity_get_databases_cli import (
    _looks_like_dna,
    _reformat_protein_fasta,
    clean_protein_fasta,
    get_databases,
)


# ---------------------------------------------------------------------------
# _looks_like_dna unit tests
# ---------------------------------------------------------------------------


class TestLooksLikeDna(unittest.TestCase):
    """Verify the DNA/RNA detection heuristic used to drop contaminated records."""

    def test_pure_dna_above_threshold_is_dna(self):
        seq = "ATGGACAACTCAATTGTTGTAGTCAGAGCTACTAAGGCC"
        self.assertTrue(_looks_like_dna(seq))

    def test_pure_rna_above_threshold_is_dna(self):
        seq = "AUGGACAACUCAAUUGUUGUAGUCAGAGCUACU"
        self.assertTrue(_looks_like_dna(seq))

    def test_dna_with_n_is_dna(self):
        seq = "ATGNNNGACAACTCAATTGTTGTAGTCAGAGCTAC"
        self.assertTrue(_looks_like_dna(seq))

    def test_lowercase_dna_is_dna(self):
        seq = "atggacaactcaattgttgtagtcagagctactaag"
        self.assertTrue(_looks_like_dna(seq))

    def test_real_protein_is_not_dna(self):
        seq = "MPKLPRGLRFGADNEILNDFQELWFPDLFIESSDTHPWYTLKGRVLNAHLDDRLPNVGGRQ"
        self.assertFalse(_looks_like_dna(seq))

    def test_short_acgt_sequence_is_kept(self):
        """Tiny ACGT peptides are biologically possible; we err on the safe
        side and only flag sequences at or above ``min_length``."""
        self.assertFalse(_looks_like_dna("ACGTACG"))

    def test_gaps_and_stops_are_ignored(self):
        seq = "ATG-GAC*AAC TCA-ATT-GTT-GTA-GTC-AGA-GCT-ACT"
        self.assertTrue(_looks_like_dna(seq))

    def test_min_length_parameter_is_respected(self):
        seq = "ATGGAC"
        self.assertFalse(_looks_like_dna(seq, min_length=20))
        self.assertTrue(_looks_like_dna(seq, min_length=5))


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------


_MIXED_FASTA = textwrap.dedent(
    """\
    >NC_139268.1
    ATGGACAACTCAATTGTTGTAGTCAGAGCTACTAAGGCCGCCTTTGTGCCAATCAAACCT
    AAATTGGAAGATGAGGTCAACTATCCTCGAGAGTTCTTTGTAGACGGAAGGATTCCTGCG
    >YP_011242554.1 polyprotein [organism=dengue virus type 2]
    MPKLPRGLRFGADNEILNDFQELWFPDLFIESSDTHPWYTLKGRVLNAHLDDRLPNVGGRQ
    IRRTPHRATVPIASSGLRPVTTVQYDPTALSFLLNARVDIRELRRELLD
    >NC_076867.1
    ATGCTATCTGCAGATGCCAGGACACGGTGGAGTGGAGCTAAACAGGACATAGAGACTCTA
    GCAAGAGGGATTAGTGGAGCCGGGAGATCAGAAGAAATCAGTTTAGATATTGAACCAGAA
    >NP_056776.2:1-3391 polyprotein [organism=zika virus]
    MAKLETVTLSNIGKDGKQTLVLNPRGVNPTNGVAALSQAGAVPALEKRVTVSVSQPSRNR
    KNYKVQVKIQNPTACTANGSCDPSVTRQAYADVTFSFTQYSTRT
    >YP_999999.1 hypothetical protein [organism=mystery virus]
    MSKTASSHNSLSAQLRRAANTRIEVEGNLALSIANDLLLAYGQSPFNSEAECISLSPRFD
    """
)


def _write_fasta(tmpdir: Path, name: str, content: str) -> Path:
    path = tmpdir / name
    path.write_text(content)
    return path


# ---------------------------------------------------------------------------
# _reformat_protein_fasta tests
# ---------------------------------------------------------------------------


class TestReformatProteinFasta(unittest.TestCase):
    """Verify DNA records are dropped during the reformat step."""

    def test_drops_dna_records_and_maps_proteins(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "protein.faa", _MIXED_FASTA)
            output_fa = tmpdir / "viral.protein.faa"

            org2taxid = {
                "dengue virus type 2": "11060",
                "zika virus": "64320",
                "mystery virus": "99999",
            }

            taxid_map = _reformat_protein_fasta(
                [input_fa],
                output_fa,
                org2taxid,
            )

            self.assertEqual(
                taxid_map,
                {
                    "YP_011242554.1": "11060",
                    "NP_056776.2": "64320",
                    "YP_999999.1": "99999",
                },
            )

            written = output_fa.read_text()
            self.assertIn(">YP_011242554.1", written)
            self.assertIn(">NP_056776.2", written)
            self.assertIn(">YP_999999.1", written)
            self.assertNotIn("NC_139268.1", written)
            self.assertNotIn("NC_076867.1", written)
            self.assertNotIn("ATGGACAACTCAATTGTTGTAGTCAGAGCTAC", written)

    def test_pure_protein_input_writes_everything(self):
        protein_only = textwrap.dedent(
            """\
            >YP_1.1 protein A [organism=virus a]
            MPKLPRGLRFGADNEILNDFQ
            >YP_2.1 protein B [organism=virus b]
            MSKTASSHNSLSAQLRRAANT
            """
        )
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "protein.faa", protein_only)
            output_fa = tmpdir / "out.faa"

            taxid_map = _reformat_protein_fasta(
                [input_fa],
                output_fa,
                {"virus a": "1", "virus b": "2"},
            )

            self.assertEqual(taxid_map, {"YP_1.1": "1", "YP_2.1": "2"})
            written = output_fa.read_text()
            self.assertIn(">YP_1.1", written)
            self.assertIn(">YP_2.1", written)


# ---------------------------------------------------------------------------
# clean-protein-fasta CLI command tests
# ---------------------------------------------------------------------------


class TestCleanProteinFastaCommand(unittest.TestCase):
    """End-to-end tests for ``viralunity get-databases clean-protein-fasta``."""

    def setUp(self):
        self.runner = CliRunner()

    def test_clean_writes_separate_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "viral.protein.faa", _MIXED_FASTA)
            output_fa = tmpdir / "viral.protein.cleaned.faa"

            result = self.runner.invoke(
                get_databases,
                [
                    "clean-protein-fasta",
                    "--input",
                    str(input_fa),
                    "--output",
                    str(output_fa),
                ],
                catch_exceptions=False,
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertTrue(output_fa.exists())

            cleaned = output_fa.read_text()
            self.assertIn(">YP_011242554.1", cleaned)
            self.assertNotIn(">NC_139268.1", cleaned)
            self.assertNotIn(">NC_076867.1", cleaned)
            self.assertIn("Kept 3 protein record(s).", result.output)
            self.assertIn("Dropped 2 DNA/RNA record(s)", result.output)

            self.assertEqual(input_fa.read_text(), _MIXED_FASTA)

    def test_clean_in_place_with_backup(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "viral.protein.faa", _MIXED_FASTA)

            result = self.runner.invoke(
                get_databases,
                [
                    "clean-protein-fasta",
                    "--input",
                    str(input_fa),
                    "--output",
                    str(input_fa),
                ],
                catch_exceptions=False,
            )

            self.assertEqual(result.exit_code, 0, result.output)

            backup = input_fa.with_suffix(input_fa.suffix + ".with_dna.bak")
            self.assertTrue(backup.exists(), "backup file should be created")
            self.assertEqual(backup.read_text(), _MIXED_FASTA)

            cleaned = input_fa.read_text()
            self.assertIn(">YP_011242554.1", cleaned)
            self.assertNotIn(">NC_139268.1", cleaned)

    def test_clean_in_place_no_backup(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "viral.protein.faa", _MIXED_FASTA)

            result = self.runner.invoke(
                get_databases,
                [
                    "clean-protein-fasta",
                    "--input",
                    str(input_fa),
                    "--output",
                    str(input_fa),
                    "--no-backup",
                ],
                catch_exceptions=False,
            )

            self.assertEqual(result.exit_code, 0, result.output)

            backup = input_fa.with_suffix(input_fa.suffix + ".with_dna.bak")
            self.assertFalse(backup.exists(), "backup must not be created with --no-backup")

            cleaned = input_fa.read_text()
            self.assertNotIn(">NC_139268.1", cleaned)
            self.assertIn(">YP_011242554.1", cleaned)

    def test_clean_protein_only_input_is_a_noop(self):
        protein_only = textwrap.dedent(
            """\
            >YP_1.1 protein A
            MPKLPRGLRFGADNEILNDFQ
            >YP_2.1 protein B
            MSKTASSHNSLSAQLRRAANT
            """
        )
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_fa = _write_fasta(tmpdir, "viral.protein.faa", protein_only)
            output_fa = tmpdir / "out.faa"

            result = self.runner.invoke(
                get_databases,
                [
                    "clean-protein-fasta",
                    "--input",
                    str(input_fa),
                    "--output",
                    str(output_fa),
                ],
                catch_exceptions=False,
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertIn("Kept 2 protein record(s).", result.output)
            self.assertIn("No DNA/RNA records found.", result.output)

    def test_clean_missing_input_fails(self):
        result = self.runner.invoke(
            get_databases,
            [
                "clean-protein-fasta",
                "--input",
                "/nonexistent/path/protein.faa",
                "--output",
                "/tmp/should_not_be_created.faa",
            ],
        )
        self.assertNotEqual(result.exit_code, 0)


if __name__ == "__main__":
    unittest.main()
