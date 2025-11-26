"""Configuration file generation for ViralUnity pipelines."""

import os
from typing import Dict, List, Any
import yaml

from viralunity.viralunity.constants import ConfigKeys, DataType
from viralunity.viralunity.exceptions import ConfigurationError


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
    
    def add_illumina_settings(
        self,
        adapters: str,
        minimum_read_length: int,
        trim: int
    ) -> None:
        """Add Illumina-specific settings to configuration.
        
        Args:
            adapters: Path to adapters file
            minimum_read_length: Minimum read length threshold
            trim: Number of bases to trim from 5' end
        """
        self.config[ConfigKeys.ADAPTERS] = adapters
        self.config[ConfigKeys.MINIMUM_LENGTH] = minimum_read_length
        self.config[ConfigKeys.TRIM] = trim
    
    def add_consensus_settings(
        self,
        reference: str,
        primer_scheme: str,
        minimum_coverage: int
    ) -> None:
        """Add consensus-specific settings to configuration.
        
        Args:
            reference: Path to reference genome
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
        os.makedirs(os.path.dirname(self.config_path), exist_ok=True)
        
        try:
            with open(self.config_path, 'w') as f:
                yaml.dump(self.config, f, default_flow_style=False, sort_keys=False)
        except (OSError, IOError) as e:
            raise ConfigurationError(
                f"Failed to write config file to {self.config_path}: {e}"
            )

