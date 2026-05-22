import tempfile
import unittest

from viralunity.config_generator import ConfigGenerator
from viralunity.constants import DataType


class TestConfigGeneratorSamples(unittest.TestCase):
    def test_nanopore_multiple_files_joined_with_space(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = f"{tmpdir}/config.yaml"
            generator = ConfigGenerator(config_path)
            paths = ["/data/a.fastq.gz", "/data/b.fastq.gz", "/data/c.fastq.gz"]
            generator.add_samples({"zika": paths}, DataType.NANOPORE)

            self.assertEqual(
                generator.config["samples"]["sample-zika"],
                "/data/a.fastq.gz /data/b.fastq.gz /data/c.fastq.gz",
            )

    def test_nanopore_single_file_unchanged(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = f"{tmpdir}/config.yaml"
            generator = ConfigGenerator(config_path)
            generator.add_samples(
                {"zika": ["/data/reads.fastq.gz"]}, DataType.NANOPORE
            )

            self.assertEqual(
                generator.config["samples"]["sample-zika"],
                "/data/reads.fastq.gz",
            )
