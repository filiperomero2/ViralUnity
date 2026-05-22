import os
import tempfile
import unittest

from viralunity.constants import DataType
from viralunity.exceptions import ValidationError
from viralunity.validators import normalize_sample_paths


class TestNormalizeSamplePaths(unittest.TestCase):
    def test_nanopore_expands_glob_in_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq_path = os.path.join(tmpdir, "reads.fastq.gz")
            with open(fastq_path, "wb"):
                pass

            samples = {"zika": [os.path.join(tmpdir, "*")]}
            normalized = normalize_sample_paths(samples, DataType.NANOPORE)

            self.assertEqual(normalized["zika"], [fastq_path])

    def test_nanopore_expands_existing_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq_path = os.path.join(tmpdir, "reads.fastq")
            with open(fastq_path, "wb"):
                pass

            samples = {"zika": [tmpdir]}
            normalized = normalize_sample_paths(samples, DataType.NANOPORE)

            self.assertEqual(normalized["zika"], [fastq_path])

    def test_nanopore_accepts_concrete_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq_path = os.path.join(tmpdir, "reads.fastq.gz")
            with open(fastq_path, "wb"):
                pass

            samples = {"zika": [fastq_path]}
            normalized = normalize_sample_paths(samples, DataType.NANOPORE)

            self.assertEqual(normalized["zika"], [fastq_path])

    def test_nanopore_accepts_multiple_matches_sorted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = []
            for name in ("b.fastq", "a.fastq"):
                path = os.path.join(tmpdir, name)
                with open(path, "wb"):
                    pass
                paths.append(path)

            samples = {"zika": [os.path.join(tmpdir, "*")]}
            normalized = normalize_sample_paths(samples, DataType.NANOPORE)

            self.assertEqual(
                normalized["zika"],
                sorted(paths),
            )

    def test_nanopore_rejects_no_matches(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            samples = {"zika": [os.path.join(tmpdir, "*")]}
            with self.assertRaises(ValidationError):
                normalize_sample_paths(samples, DataType.NANOPORE)

    def test_illumina_requires_two_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            r1 = os.path.join(tmpdir, "r1.fastq.gz")
            r2 = os.path.join(tmpdir, "r2.fastq.gz")
            for path in (r1, r2):
                with open(path, "wb"):
                    pass

            samples = {"sample1": [r1, r2]}
            normalized = normalize_sample_paths(samples, DataType.ILLUMINA)

            self.assertEqual(normalized["sample1"], [r1, r2])
