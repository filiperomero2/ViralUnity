"""Configuration file generation for ViralUnity pipelines."""

import os
from typing import Dict, List, Any, Optional
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
    
    def add_illumina_settings(
        self,
        adapters: str,
        minimum_read_length: int,
        trim: int,
        trim_head: Optional[int] = None,
        trim_tail: Optional[int] = None,
        cut_front_mean_quality: int = 20,
        cut_tail_mean_quality: int = 20,
        cut_right_window_size: int = 4,
        cut_right_mean_quality: int = 20,
    ) -> None:
        """Add Illumina-specific settings to configuration (fastp QC).

        Args:
            adapters: Path to adapters file or "NA" for auto-detection
            minimum_read_length: Minimum read length threshold
            trim: Number of bases to trim from 5' end (used as trim_head if trim_head not set)
            trim_head: Bases to trim from 5' (overrides trim if set)
            trim_tail: Bases to trim from 3'
            cut_front_mean_quality: fastp cut_front mean quality threshold
            cut_tail_mean_quality: fastp cut_tail mean quality threshold
            cut_right_window_size: fastp cut_right window size
            cut_right_mean_quality: fastp cut_right mean quality threshold
        """
        self.config[ConfigKeys.ADAPTERS] = adapters
        self.config[ConfigKeys.MINIMUM_LENGTH] = minimum_read_length
        self.config[ConfigKeys.TRIM] = trim
        head = trim_head if trim_head is not None else trim
        self.config[ConfigKeys.TRIM_HEAD] = head
        self.config[ConfigKeys.TRIM_TAIL] = trim_tail if trim_tail is not None else 0
        self.config[ConfigKeys.CUT_FRONT_MEAN_QUALITY] = cut_front_mean_quality
        self.config[ConfigKeys.CUT_TAIL_MEAN_QUALITY] = cut_tail_mean_quality
        self.config[ConfigKeys.CUT_RIGHT_WINDOW_SIZE] = cut_right_window_size
        self.config[ConfigKeys.CUT_RIGHT_MEAN_QUALITY] = cut_right_mean_quality
    
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
        remove_unclassified_reads: bool,
        host_reference: str = "NA",
        deacon_index: str = "NA",
        taxdump: str = "NA",
        run_denovo_assembly: bool = False,
        run_kraken2_reads: bool = True,
        run_kraken2_contigs: bool = True,
        run_diamond_reads: bool = False,
        run_diamond_contigs: bool = False,
        assembly_summary: str = "NA",
        diamond_database: str = "NA",
        diamond_sensitivity: str = "sensitive",
        evalue: float = 0.001,
        bleed_fraction: float = 0.005,
        negative_controls: Optional[List[str]] = None,
        negative_p_threshold: float = 0.01,
        minimum_hit_group: int = 4,
    ) -> None:
        """Add metagenomics-specific settings to configuration.

        Args:
            kraken2_database: Path to Kraken2 database
            krona_database: Path to Krona database
            remove_human_reads: Whether to remove human reads
            remove_unclassified_reads: Whether to remove unclassified reads
            host_reference: Path to host genome FASTA for dehosting (or "NA")
            deacon_index: Path to Deacon minimizer index for host depletion (or "NA"). If set, used instead of host_reference for dehosting.
            taxdump: Path to NCBI taxdump dir (nodes.dmp, names.dmp)
            run_denovo_assembly: Whether to run MEGAHIT assembly
            run_kraken2_reads: Whether to run Kraken2 on reads
            run_kraken2_contigs: Whether to run Kraken2 on contigs (if assembly)
            run_diamond_reads: Whether to run DIAMOND on reads
            run_diamond_contigs: Whether to run DIAMOND on contigs (if assembly)
            assembly_summary: NCBI assembly summary for Diamond taxonomy (or "NA")
            diamond_database: Path to Diamond DB (protein FASTA)
            diamond_sensitivity: Diamond sensitivity (e.g. sensitive, mid-sensitive)
            evalue: E-value threshold for Diamond
            bleed_fraction: Max-RPM bleed filter fraction
            negative_controls: Sample IDs to use as negative controls
            negative_p_threshold: p-value threshold for negative filter
            minimum_hit_group: Kraken2 --minimum-hit-group (default: 4)
        """
        self.config[ConfigKeys.KRAKEN2_DATABASE] = kraken2_database
        self.config[ConfigKeys.KRONA_DATABASE] = krona_database
        self.config[ConfigKeys.REMOVE_HUMAN_READS] = remove_human_reads
        self.config[ConfigKeys.REMOVE_UNCLASSIFIED_READS] = remove_unclassified_reads
        self.config[ConfigKeys.HOST_REFERENCE] = host_reference
        self.config[ConfigKeys.DEACON_INDEX] = deacon_index
        self.config[ConfigKeys.TAXDUMP] = taxdump
        self.config[ConfigKeys.RUN_DENOVO_ASSEMBLY] = run_denovo_assembly
        self.config[ConfigKeys.RUN_KRAKEN2_READS] = run_kraken2_reads
        self.config[ConfigKeys.RUN_KRAKEN2_CONTIGS] = run_kraken2_contigs
        self.config[ConfigKeys.RUN_DIAMOND_READS] = run_diamond_reads
        self.config[ConfigKeys.RUN_DIAMOND_CONTIGS] = run_diamond_contigs
        self.config[ConfigKeys.ASSEMBLY_SUMMARY] = assembly_summary
        self.config[ConfigKeys.DIAMOND_DATABASE] = diamond_database
        self.config[ConfigKeys.DIAMOND_SENSITIVITY] = diamond_sensitivity
        self.config[ConfigKeys.EVALUE] = evalue
        self.config[ConfigKeys.BLEED_FRACTION] = bleed_fraction
        self.config[ConfigKeys.NEGATIVE_CONTROLS] = negative_controls or []
        self.config[ConfigKeys.NEGATIVE_P_THRESHOLD] = negative_p_threshold
        self.config[ConfigKeys.MINIMUM_HIT_GROUP] = minimum_hit_group
    
    def add_nanopore_settings(
        self,
        run_polish_racon: bool = False,
        run_polish_medaka: bool = False,
        medaka_model: Optional[str] = None,
    ) -> None:
        """Add Nanopore-specific settings (polishing: Racon, Medaka).
        
        Args:
            run_polish_racon: Whether to run Racon polishing after MEGAHIT.
            run_polish_medaka: Whether to run Medaka polishing (after Racon if both enabled).
            medaka_model: Medaka model name (e.g. r941_min_high_g360). Optional; Medaka uses default if not set.
        """
        self.config[ConfigKeys.RUN_POLISH_RACON] = run_polish_racon
        self.config[ConfigKeys.RUN_POLISH_MEDAKA] = run_polish_medaka
        if medaka_model is not None:
            self.config[ConfigKeys.MEDAKA_MODEL] = medaka_model
    
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

