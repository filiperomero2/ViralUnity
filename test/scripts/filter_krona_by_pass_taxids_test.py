"""Tests for viralunity.scripts.python.filter_krona_by_pass_taxids."""

import os
import tempfile
import unittest

import pandas as pd

from viralunity.scripts.python.filter_krona_by_pass_taxids import (
    _lineage_passes,
    build_pass_taxids,
    filter_krona_input,
    get_lineage,
    load_taxdump,
    run,
)

# Synthetic taxonomy used by all tests:
#     1 (no rank, root)
#       └── 1000 (family)
#             └── 2000 (genus)
#                   ├── 3001 (species)
#                   │     └── 4001 (no rank / strain)
#                   └── 3002 (species)
#                         └── 4002 (no rank / strain)
#       └── 1100 (family)
#             └── 2100 (genus)
#                   └── 3100 (species)

NODES_RECORDS = [
    ("1", "1", "no rank"),
    ("1000", "1", "family"),
    ("2000", "1000", "genus"),
    ("3001", "2000", "species"),
    ("3002", "2000", "species"),
    ("4001", "3001", "no rank"),
    ("4002", "3002", "no rank"),
    ("1100", "1", "family"),
    ("2100", "1100", "genus"),
    ("3100", "2100", "species"),
]

NAMES_RECORDS = [
    ("1", "root", "", "scientific name"),
    ("1000", "FamilyA", "", "scientific name"),
    ("2000", "GenusA", "", "scientific name"),
    ("3001", "SpeciesA1", "", "scientific name"),
    ("3002", "SpeciesA2", "", "scientific name"),
    ("4001", "StrainA1a", "", "scientific name"),
    ("4002", "StrainA2a", "", "scientific name"),
    ("1100", "FamilyB", "", "scientific name"),
    ("2100", "GenusB", "", "scientific name"),
    ("3100", "SpeciesB1", "", "scientific name"),
]


def _write_taxdump(directory: str) -> str:
    nodes_path = os.path.join(directory, "nodes.dmp")
    names_path = os.path.join(directory, "names.dmp")
    with open(nodes_path, "w") as f:
        for taxid, parent, rank in NODES_RECORDS:
            f.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t-\t|\n")
    with open(names_path, "w") as f:
        for taxid, name, _unique, name_class in NAMES_RECORDS:
            f.write(f"{taxid}\t|\t{name}\t|\t\t|\t{name_class}\t|\n")
    return directory


def _write_krona_input(path: str, rows):
    with open(path, "w") as f:
        for query, taxid in rows:
            f.write(f"{query}\t{taxid}\n")


def _read_lines(path: str):
    with open(path) as f:
        return [line.rstrip("\n") for line in f if line.strip()]


def _summary_rows():
    """A two-sample, two-context summary table with some passing taxa."""
    return pd.DataFrame(
        [
            # sample A, kraken2/reads: SpeciesA1 passes
            {
                "sample": "A",
                "tool": "kraken2",
                "mode": "reads",
                "rank": "species",
                "taxid": "3001",
                "bleed_pass": True,
            },
            # sample A, kraken2/reads: SpeciesA2 fails
            {
                "sample": "A",
                "tool": "kraken2",
                "mode": "reads",
                "rank": "species",
                "taxid": "3002",
                "bleed_pass": False,
            },
            # sample A, kraken2/reads: GenusB passes (genus-rank entry)
            {
                "sample": "A",
                "tool": "kraken2",
                "mode": "reads",
                "rank": "genus",
                "taxid": "2100",
                "bleed_pass": True,
            },
            # sample B, kraken2/reads: SpeciesB1 passes — must NOT affect sample A
            {
                "sample": "B",
                "tool": "kraken2",
                "mode": "reads",
                "rank": "species",
                "taxid": "3100",
                "bleed_pass": True,
            },
            # sample A, kraken2/contigs: nothing passes (different mode)
            {
                "sample": "A",
                "tool": "kraken2",
                "mode": "contigs",
                "rank": "species",
                "taxid": "3001",
                "bleed_pass": False,
            },
            # sample A, diamond/reads: SpeciesA2 passes — must NOT affect kraken2/reads
            {
                "sample": "A",
                "tool": "diamond",
                "mode": "reads",
                "rank": "species",
                "taxid": "3002",
                "bleed_pass": True,
            },
        ]
    )


class TestLineageHelpers(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        _write_taxdump(self._tmp.name)
        self.parent_map, self.rank_map = load_taxdump(
            os.path.join(self._tmp.name, "nodes.dmp"),
            os.path.join(self._tmp.name, "names.dmp"),
        )

    def test_load_taxdump_parses_parents_and_ranks(self):
        self.assertEqual(self.parent_map["4001"], "3001")
        self.assertEqual(self.parent_map["3001"], "2000")
        self.assertEqual(self.rank_map["3001"], "species")
        self.assertEqual(self.rank_map["2000"], "genus")
        self.assertEqual(self.rank_map["1000"], "family")

    def test_get_lineage_includes_self_and_root(self):
        lineage = get_lineage("4001", self.parent_map)
        self.assertEqual(lineage, ["4001", "3001", "2000", "1000", "1"])

    def test_lineage_passes_when_species_in_pass_set(self):
        self.assertTrue(_lineage_passes("4001", {"3001"}, self.parent_map, self.rank_map))

    def test_lineage_passes_when_genus_in_pass_set(self):
        # Strain 4002 lives under SpeciesA2 / GenusA / FamilyA. If only genus
        # 2000 passes, the strain row must still be kept.
        self.assertTrue(_lineage_passes("4002", {"2000"}, self.parent_map, self.rank_map))

    def test_lineage_fails_when_no_relevant_rank_passes(self):
        self.assertFalse(_lineage_passes("4001", {"3002"}, self.parent_map, self.rank_map))

    def test_taxid_zero_never_passes(self):
        self.assertFalse(_lineage_passes("0", {"3001", "3002"}, self.parent_map, self.rank_map))

    def test_unknown_taxid_never_passes(self):
        self.assertFalse(_lineage_passes("999999", {"3001"}, self.parent_map, self.rank_map))


class TestBuildPassTaxids(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.summary_path = os.path.join(self._tmp.name, "summary.tsv")
        _summary_rows().to_csv(self.summary_path, sep="\t", index=False)

    def test_pass_set_is_restricted_by_sample(self):
        a = build_pass_taxids(self.summary_path, sample="A", tool="kraken2", mode="reads")
        b = build_pass_taxids(self.summary_path, sample="B", tool="kraken2", mode="reads")
        self.assertEqual(a, {"3001", "2100"})
        self.assertEqual(b, {"3100"})

    def test_pass_set_is_restricted_by_mode(self):
        passes = build_pass_taxids(self.summary_path, sample="A", tool="kraken2", mode="contigs")
        self.assertEqual(passes, set())

    def test_pass_set_is_restricted_by_tool(self):
        kraken = build_pass_taxids(self.summary_path, sample="A", tool="kraken2", mode="reads")
        diamond = build_pass_taxids(self.summary_path, sample="A", tool="diamond", mode="reads")
        self.assertEqual(kraken, {"3001", "2100"})
        self.assertEqual(diamond, {"3002"})

    def test_neg_pass_drops_false_keeps_true_and_na(self):
        df = _summary_rows()
        # bleed_pass is True for these rows; vary neg_pass.
        df["neg_pass"] = [True, False, pd.NA, True, False, False]
        df.to_csv(self.summary_path, sep="\t", index=False)
        # For sample A, kraken2/reads we have: row 0 (3001, neg=True keep),
        # row 1 (3002, bleed_pass=False already drops), row 2 (2100, neg=NA keep).
        passes = build_pass_taxids(self.summary_path, sample="A", tool="kraken2", mode="reads")
        self.assertEqual(passes, {"3001", "2100"})

    def test_missing_taxid_column_raises(self):
        bad = pd.DataFrame(
            [{"sample": "A", "tool": "kraken2", "mode": "reads", "bleed_pass": True}]
        )
        bad.to_csv(self.summary_path, sep="\t", index=False)
        with self.assertRaises(ValueError):
            build_pass_taxids(self.summary_path, sample="A", tool="kraken2", mode="reads")

    def test_missing_bleed_pass_column_raises(self):
        bad = pd.DataFrame([{"sample": "A", "tool": "kraken2", "mode": "reads", "taxid": "3001"}])
        bad.to_csv(self.summary_path, sep="\t", index=False)
        with self.assertRaises(ValueError):
            build_pass_taxids(self.summary_path, sample="A", tool="kraken2", mode="reads")


class TestFilterKronaInput(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        _write_taxdump(self._tmp.name)
        self.parent_map, self.rank_map = load_taxdump(
            os.path.join(self._tmp.name, "nodes.dmp"),
            os.path.join(self._tmp.name, "names.dmp"),
        )
        self.krona_in = os.path.join(self._tmp.name, "krona.tsv")
        self.krona_out = os.path.join(self._tmp.name, "filtered.tsv")

    def test_keeps_strain_when_species_passes(self):
        _write_krona_input(
            self.krona_in,
            [
                ("contig_strain", "4001"),
                ("contig_species_a2", "3002"),
                ("contig_unknown", "9999"),
                ("contig_zero", "0"),
            ],
        )
        kept, total = filter_krona_input(
            self.krona_in,
            self.krona_out,
            {"3001"},
            self.parent_map,
            self.rank_map,
        )
        self.assertEqual(total, 4)
        self.assertEqual(kept, 1)
        self.assertEqual(_read_lines(self.krona_out), ["contig_strain\t4001"])

    def test_keeps_species_direct_match(self):
        _write_krona_input(
            self.krona_in,
            [
                ("contig_species_a1", "3001"),
                ("contig_species_a2", "3002"),
            ],
        )
        kept, _ = filter_krona_input(
            self.krona_in,
            self.krona_out,
            {"3001"},
            self.parent_map,
            self.rank_map,
        )
        self.assertEqual(kept, 1)
        self.assertEqual(_read_lines(self.krona_out), ["contig_species_a1\t3001"])

    def test_keeps_when_genus_passes(self):
        _write_krona_input(
            self.krona_in,
            [
                ("contig_strain_a1a", "4001"),
                ("contig_strain_a2a", "4002"),
                ("contig_species_b1", "3100"),
            ],
        )
        kept, _ = filter_krona_input(
            self.krona_in,
            self.krona_out,
            {"2000"},
            self.parent_map,
            self.rank_map,
        )
        # 4001 and 4002 are both under genus 2000; 3100 is under genus 2100.
        self.assertEqual(kept, 2)
        self.assertEqual(
            sorted(_read_lines(self.krona_out)),
            sorted(["contig_strain_a1a\t4001", "contig_strain_a2a\t4002"]),
        )

    def test_empty_pass_set_yields_empty_output(self):
        _write_krona_input(self.krona_in, [("c1", "3001"), ("c2", "3002")])
        kept, _ = filter_krona_input(
            self.krona_in,
            self.krona_out,
            set(),
            self.parent_map,
            self.rank_map,
        )
        self.assertEqual(kept, 0)
        self.assertEqual(_read_lines(self.krona_out), [])

    def test_empty_krona_input_yields_empty_output(self):
        _write_krona_input(self.krona_in, [])
        kept, total = filter_krona_input(
            self.krona_in,
            self.krona_out,
            {"3001"},
            self.parent_map,
            self.rank_map,
        )
        self.assertEqual((kept, total), (0, 0))
        self.assertEqual(_read_lines(self.krona_out), [])


class TestRunEndToEnd(unittest.TestCase):
    """`run()` glues build_pass_taxids + load_taxdump + filter_krona_input."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        _write_taxdump(self._tmp.name)
        self.summary = os.path.join(self._tmp.name, "summary.tsv")
        _summary_rows().to_csv(self.summary, sep="\t", index=False)
        self.krona_in = os.path.join(self._tmp.name, "krona.tsv")
        self.krona_out = os.path.join(self._tmp.name, "filtered.tsv")
        self.nodes = os.path.join(self._tmp.name, "nodes.dmp")
        self.names = os.path.join(self._tmp.name, "names.dmp")

    def test_run_lineage_aware_filtering(self):
        # Sample A, kraken2/reads pass set is {3001, 2100}. Therefore:
        #   contig under species 3001 → kept (species match)
        #   contig under strain 4001 (descendant of species 3001) → kept (lineage)
        #   contig under species 3002 → dropped
        #   contig under species 3100 (descendant of genus 2100) → kept (genus match)
        _write_krona_input(
            self.krona_in,
            [
                ("c_3001", "3001"),
                ("c_4001", "4001"),
                ("c_3002", "3002"),
                ("c_3100", "3100"),
            ],
        )
        run(
            summary=self.summary,
            krona_input=self.krona_in,
            output=self.krona_out,
            sample="A",
            tool="kraken2",
            mode="reads",
            nodes_dmp=self.nodes,
            names_dmp=self.names,
        )
        self.assertEqual(
            sorted(_read_lines(self.krona_out)),
            sorted(["c_3001\t3001", "c_4001\t4001", "c_3100\t3100"]),
        )

    def test_run_sample_isolation(self):
        # Sample B's pass set ({3100}) must not pull in sample A's hits.
        _write_krona_input(
            self.krona_in,
            [
                ("c_3001", "3001"),
                ("c_3100", "3100"),
            ],
        )
        run(
            summary=self.summary,
            krona_input=self.krona_in,
            output=self.krona_out,
            sample="B",
            tool="kraken2",
            mode="reads",
            nodes_dmp=self.nodes,
            names_dmp=self.names,
        )
        self.assertEqual(_read_lines(self.krona_out), ["c_3100\t3100"])

    def test_run_tool_mode_isolation(self):
        # Sample A, diamond/reads pass set is {3002}; the kraken2/reads pass
        # set must not leak across.
        _write_krona_input(
            self.krona_in,
            [
                ("c_3001", "3001"),
                ("c_3002", "3002"),
            ],
        )
        run(
            summary=self.summary,
            krona_input=self.krona_in,
            output=self.krona_out,
            sample="A",
            tool="diamond",
            mode="reads",
            nodes_dmp=self.nodes,
            names_dmp=self.names,
        )
        self.assertEqual(_read_lines(self.krona_out), ["c_3002\t3002"])

    def test_run_with_neg_pass_false_drops_taxon(self):
        df = _summary_rows()
        df["neg_pass"] = [False, False, True, True, False, False]
        df.to_csv(self.summary, sep="\t", index=False)
        # For sample A, kraken2/reads pass set with bleed+neg becomes {2100}.
        # contig under species 3001 (no longer passing) must be dropped;
        # contig under species 3100 (under passing genus 2100) must be kept.
        _write_krona_input(
            self.krona_in,
            [
                ("c_3001", "3001"),
                ("c_3100", "3100"),
            ],
        )
        run(
            summary=self.summary,
            krona_input=self.krona_in,
            output=self.krona_out,
            sample="A",
            tool="kraken2",
            mode="reads",
            nodes_dmp=self.nodes,
            names_dmp=self.names,
        )
        self.assertEqual(_read_lines(self.krona_out), ["c_3100\t3100"])


if __name__ == "__main__":
    unittest.main()
