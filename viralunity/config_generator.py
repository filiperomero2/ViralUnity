"""Configuration file generation for ViralUnity pipelines."""

import os
from typing import Dict, List, Any, Union
import yaml

from viralunity.constants import ConfigKeys, DataType
from viralunity.exceptions import ConfigurationError


class ConfigGenerator:
    """Generates YAML configuration files for Snakemake workflows."""
    
    def __init__(self, config_path: str):
        """Initialize config generator.
        
        Args:
            config_path: Path where the config file will be written
        """
        self.config_path = config_path
        self.config: Dict[str, Any] = {}
    
    def add_samples(self, samples: Dict[str, List[str]], data_type: str) -> None:
        """Add samples to configuration.
        
        Args:
            samples: Dictionary mapping sample names to file paths
            data_type: Type of sequencing data (illumina or nanopore)
        """
        formatted_samples = {}
        for sample_name, file_paths in samples.items():
            key = f"sample-{sample_name}"
            if data_type == DataType.ILLUMINA:
                if len(file_paths) != 2:
                    raise ConfigurationError(
                        f"Illumina sample {sample_name} must have 2 files, "
                        f"found {len(file_paths)}"
                    )
                formatted_samples[key] = f"{file_paths[0]} {file_paths[1]}"
            else:
                if len(file_paths) != 1:
                    raise ConfigurationError(
                        f"Nanopore sample {sample_name} must have 1 file, "
                        f"found {len(file_paths)}"
                    )
                formatted_samples[key] = file_paths[0]
        
        self.config[ConfigKeys.SAMPLES] = formatted_samples
        self.config[ConfigKeys.DATA] = data_type
    
    def add_output(self, output_dir: str, run_name: str) -> None:
        """Add output directory to configuration.
        
        Args:
            output_dir: Base output directory
            run_name: Name of the run
        """
        self.config[ConfigKeys.OUTPUT] = os.path.join(output_dir, run_name, "")
    
    def add_threads(self, threads: int) -> None:
        """Add thread count to configuration.
        
        Args:
            threads: Number of threads
        """
        self.config[ConfigKeys.THREADS] = threads

    def add_nanopore_settings(
        self,
        minimum_read_length: int,
        af_threshold: float,
        chunk_size: int,
        clair3_model: str,
        variant_quality: int,
        variant_depth: int,
        minimum_map_quality: int
    ) -> None:
        """Add Nanopore-specific settings to configuration.

        Args:
            minimum_read_length: Minimum read length threshold
            af_threshold: Allele frequency threshold to call a variant into consensus
            chunk_size: Size of chunks to process [clair3]
            clair3_model: Model to use for variant calling [clair3]
            variant_quality: Minimum variant quality to call a variant into consensus [clair3]
            variant_depth: Minimum alt allele depth to call a variant into consensus [clair3]
            minimum_map_quality: Minimum map quality to call a variant into consensus [clair3]
        """
        self.config[ConfigKeys.MINIMUM_LENGTH] = minimum_read_length
        self.config[ConfigKeys.AF_THRESHOLD] = af_threshold
        self.config[ConfigKeys.CHUNK_SIZE] = chunk_size
        self.config[ConfigKeys.CLAIR3_MODEL] = clair3_model
        self.config[ConfigKeys.VARIANT_QUALITY] = variant_quality
        self.config[ConfigKeys.VARIANT_DEPTH] = variant_depth
        self.config[ConfigKeys.MINIMUM_MAP_QUALITY] = minimum_map_quality

    def add_illumina_settings(
        self,
        adapters: str,
        minimum_read_length: int,
        trim_head: int,
        trim_tail: int,
        cut_front_mean_quality: int = 20,
        cut_tail_mean_quality: int = 20,
        cut_right_window_size: int = 4,
        cut_right_mean_quality: int = 20,
        af_threshold: float = 0.5,
        af_isnv_threshold: float = 0.05,
        run_isnv: bool = False,
    ) -> None:
        """Add Illumina-specific settings to configuration.
        
        Args:
            adapters: Path to adapters file
            minimum_read_length: Minimum read length threshold
            trim_head: Number of bases to trim from 5' end
            trim_tail: Number of bases to trim from 3' end
            cut_front_mean_quality: Mean quality requirement option for cut_front
            cut_tail_mean_quality: Mean quality requirement option for cut_tail
            cut_right_window_size: Window size for cut_right
            cut_right_mean_quality: Mean quality requirement option for cut_right
            af_threshold: Allele frequency threshold to call a variant into consensus
            af_isnv_threshold: Minimum allele frequency threshold to call a variant into iSNV analysis
        """
        self.config[ConfigKeys.ADAPTERS] = adapters
        self.config[ConfigKeys.MINIMUM_LENGTH] = minimum_read_length
        self.config[ConfigKeys.TRIM_HEAD] = trim_head
        self.config[ConfigKeys.TRIM_TAIL] = trim_tail
        self.config[ConfigKeys.CUT_FRONT_MEAN_QUALITY] = cut_front_mean_quality
        self.config[ConfigKeys.CUT_TAIL_MEAN_QUALITY] = cut_tail_mean_quality
        self.config[ConfigKeys.CUT_RIGHT_WINDOW_SIZE] = cut_right_window_size
        self.config[ConfigKeys.CUT_RIGHT_MEAN_QUALITY] = cut_right_mean_quality    
        self.config[ConfigKeys.AF_THRESHOLD] = af_threshold
        self.config[ConfigKeys.AF_ISNV_THRESHOLD] = af_isnv_threshold
        self.config[ConfigKeys.RUN_ISNV] = run_isnv
    
    def add_consensus_settings(
        self,
        reference: Union[str, Dict[str, str]],
        primer_scheme: str,
        minimum_coverage: int
    ) -> None:
        """Add consensus-specific settings to configuration.
        
        Args:
            reference: Path to reference genome (str) or dict mapping
                segment names to paths for segmented viruses
            primer_scheme: Path to primer scheme file or "NA"
            minimum_coverage: Minimum coverage for consensus
        """
        self.config[ConfigKeys.REFERENCE] = reference
        self.config[ConfigKeys.SCHEME] = primer_scheme
        self.config[ConfigKeys.MINIMUM_DEPTH] = minimum_coverage
    
    def add_metagenomics_settings(
        self,
        kraken2_database: str,
        krona_database: str,
        remove_human_reads: bool,
        remove_unclassified_reads: bool
    ) -> None:
        """Add metagenomics-specific settings to configuration.
        
        Args:
            kraken2_database: Path to Kraken2 database
            krona_database: Path to Krona database
            remove_human_reads: Whether to remove human reads
            remove_unclassified_reads: Whether to remove unclassified reads
        """
        self.config[ConfigKeys.KRAKEN2_DATABASE] = kraken2_database
        self.config[ConfigKeys.KRONA_DATABASE] = krona_database
        self.config[ConfigKeys.REMOVE_HUMAN_READS] = remove_human_reads
        self.config[ConfigKeys.REMOVE_UNCLASSIFIED_READS] = remove_unclassified_reads
    
    def add_workflow_path(self, workflow_path: str) -> None:
        """Add workflow path to configuration.
        
        Args:
            workflow_path: Path to the workflow directory
        """
        self.config["workflow_path"] = workflow_path
    
    def save(self) -> None:
        """Save configuration to YAML file.
        
        Raises:
            ConfigurationError: If config directory cannot be created
        """
        config_dir = os.path.dirname(self.config_path)
        if config_dir:  # Only create directory if path contains a directory component
            os.makedirs(config_dir, exist_ok=True)
        
        try:
            with open(self.config_path, 'w') as f:
                yaml.dump(self.config, f, default_flow_style=False, sort_keys=False)
        except (OSError, IOError) as e:
            raise ConfigurationError(
                f"Failed to write config file to {self.config_path}: {e}"
            )

