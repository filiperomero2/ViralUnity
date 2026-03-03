"""Validation functions for ViralUnity pipeline arguments and data."""

import os
from typing import Dict, Any, List
import pandas as pd

from viralunity.exceptions import (
    Kraken2DatabaseNotFoundError,
    KronaDatabaseNotFoundError,
    ValidationError,
    FileNotFoundError,
    SampleSheetError,
    SampleConfigurationNotFoundError,
    AdaptersNotFoundError,
    ReferenceNotFoundError,
    PrimerSchemeNotFoundError,
    TaxdumpNotFoundError,
    DiamondDatabaseNotFoundError,
)
from viralunity.constants import DataType


def validate_file_exists(file_path: str, description: str = "File") -> None:
    """Validate that a file exists.
    
    Args:
        file_path: Path to the file
        description: Description of the file for error messages
        
    Raises:
        FileNotFoundError: If the file does not exist
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"{description} does not exist: {file_path}")


def validate_directory_exists(directory_path: str, description: str = "Directory") -> None:
    """Validate that a directory exists.
    
    Args:
        directory_path: Path to the directory
        description: Description of the directory for error messages
        
    Raises:
        FileNotFoundError: If the directory does not exist
    """
    if not os.path.isdir(directory_path):
        raise FileNotFoundError(f"{description} does not exist: {directory_path}")


def validate_sample_sheet(
    sample_sheet_path: str,
    data_type: str
) -> Dict[str, List[str]]:
    """Validate and parse sample sheet file.
    
    Args:
        sample_sheet_path: Path to the sample sheet CSV file
        data_type: Type of sequencing data (illumina or nanopore)
        
    Returns:
        Dictionary mapping sample names to file paths
        
    Raises:
        SampleSheetError: If the sample sheet is invalid
        FileNotFoundError: If the sample sheet file does not exist
    """
    validate_file_exists(sample_sheet_path, "Sample sheet file")
    
    try:
        df = pd.read_csv(sample_sheet_path, header=None)
    except Exception as e:
        raise SampleSheetError(f"Failed to read sample sheet: {e}")
    
    samples = {}
    
    for idx in df.index:
        sample_name = df.iloc[idx, 0]
        
        if data_type == DataType.ILLUMINA:
            if df.shape[1] < 3:
                raise SampleSheetError(
                    f"Illumina sample sheet must have at least 3 columns. "
                    f"Found {df.shape[1]} columns."
                )
            
            sample_r1 = df.iloc[idx, 1]
            sample_r2 = df.iloc[idx, 2]
            
            validate_file_exists(sample_r1, f"R1 file for sample {sample_name}")
            validate_file_exists(sample_r2, f"R2 file for sample {sample_name}")
            
            samples[sample_name] = [sample_r1, sample_r2]
            
        else:  # nanopore
            if df.shape[1] < 2:
                raise SampleSheetError(
                    f"Nanopore sample sheet must have at least 2 columns. "
                    f"Found {df.shape[1]} columns."
                )
            
            sample_file = df.iloc[idx, 1]
            validate_file_exists(sample_file, f"File for sample {sample_name}")
            
            samples[sample_name] = [sample_file]
    
    if not samples:
        raise SampleSheetError("No valid samples found in sample sheet")
    
    return samples


def validate_illumina_requirements(args: Dict[str, Any]) -> None:
    """Validate Illumina-specific requirements (fastp: adapters optional)."""
    if args.get("data_type") != DataType.ILLUMINA:
        return

    adapters = args.get("adapters")
    if adapters and str(adapters).strip() and str(adapters).strip() != "NA":
        try:
            validate_file_exists(adapters, "Illumina adapter sequences file")
        except FileNotFoundError as e:
            raise AdaptersNotFoundError(f"Illumina adapter sequences file does not exist: {e}")


def validate_nanopore_requirements(args: Dict[str, Any]) -> None:
    """Validate Nanopore-specific requirements (e.g. polishing options)."""
    if args.get("data_type") != DataType.NANOPORE:
        return
    # Optional: validate medaka_model if provided (Medaka accepts known model names)
    # For now no strict validation; extend as needed.


def validate_consensus_requirements(args: Dict[str, Any]) -> None:
    """Validate consensus pipeline requirements.
    
    Args:
        args: Dictionary of pipeline arguments
        
    Raises:
        ValidationError: If consensus requirements are not met
    """
    reference = args.get("reference")
    if not reference:
        raise ValidationError("Reference sequence file is required")
    
    try:
        validate_file_exists(reference, "Reference sequence file")
    except FileNotFoundError as e:
        raise ReferenceNotFoundError(f"Reference sequence file does not exist: {e}")
    
    primer_scheme = args.get("primer_scheme")
    if primer_scheme:
        try:
            validate_file_exists(primer_scheme, "Primer scheme file")
        except FileNotFoundError as e:
            raise PrimerSchemeNotFoundError(f"Primer scheme file does not exist: {e}")


def validate_metagenomics_requirements(args: Dict[str, Any]) -> None:
    """Validate metagenomics pipeline requirements.

    Kraken2 and Diamond are optional; validate only the resources for the tools
    the user has enabled (run_kraken2_reads, run_kraken2_contigs, run_diamond_reads,
    run_diamond_contigs).
    """
    run_k2_reads = args.get("run_kraken2_reads", True)
    run_k2_contigs = args.get("run_kraken2_contigs", True)
    run_diamond_reads = args.get("run_diamond_reads", False)
    run_diamond_contigs = args.get("run_diamond_contigs", False)
    any_kraken2 = run_k2_reads or run_k2_contigs
    any_diamond = run_diamond_reads or run_diamond_contigs
    any_classification = any_kraken2 or any_diamond

    # Krona database: required whenever any classification is run (Krona plots for both tools)
    if any_classification:
        krona_db = args.get("krona_database")
        if not krona_db or krona_db == "NA":
            raise KronaDatabaseNotFoundError(
                "Krona database directory is required when running any classification "
                "(Kraken2 and/or Diamond). Set --krona-database."
            )
        try:
            validate_directory_exists(krona_db, "Krona database directory")
        except FileNotFoundError as e:
            raise KronaDatabaseNotFoundError(f"Krona database directory does not exist: {e}")

    # Kraken2 database: required only when Kraken2 is enabled
    if any_kraken2:
        kraken2_db = args.get("kraken2_database")
        if not kraken2_db or kraken2_db == "NA":
            raise Kraken2DatabaseNotFoundError(
                "Kraken2 database directory is required when running Kraken2. "
                "Set --kraken2-database or disable Kraken2 with --no-kraken2-reads / --no-kraken2-contigs."
            )
        try:
            validate_directory_exists(kraken2_db, "Kraken2 database directory")
        except FileNotFoundError as e:
            raise Kraken2DatabaseNotFoundError(f"Kraken2 database directory does not exist: {e}")

    # Taxdump: required for taxonomic summaries whenever any classification is run
    if any_classification:
        taxdump = args.get("taxdump", "NA")
        if not taxdump or taxdump == "NA":
            raise TaxdumpNotFoundError(
                "taxdump directory (NCBI nodes.dmp, names.dmp) is required for taxonomic summaries "
                "when running any classification. Set --taxdump."
            )
        try:
            validate_directory_exists(taxdump, "Taxdump directory")
        except FileNotFoundError as e:
            raise TaxdumpNotFoundError(f"Taxdump directory does not exist: {e}")
        nodes = os.path.join(taxdump, "nodes.dmp")
        names = os.path.join(taxdump, "names.dmp")
        if not os.path.isfile(nodes) or not os.path.isfile(names):
            raise TaxdumpNotFoundError(
                f"Taxdump directory must contain nodes.dmp and names.dmp: {taxdump}"
            )

    # Diamond: require database and assembly summary only when Diamond is enabled
    if any_diamond:
        diamond_db = args.get("diamond_database", "NA")
        if not diamond_db or diamond_db == "NA":
            raise DiamondDatabaseNotFoundError(
                "diamond_database is required when running Diamond. "
                "Set --diamond-database or do not use --run-diamond-reads / --run-diamond-contigs."
            )
        assembly = args.get("assembly_summary", "NA")
        if not assembly or assembly == "NA":
            raise DiamondDatabaseNotFoundError(
                "assembly_summary is required when running Diamond (for taxonomy mapping). Set --assembly-summary."
            )
        if not os.path.isfile(assembly) and not (assembly.endswith(".gz") and os.path.isfile(assembly)):
            raise DiamondDatabaseNotFoundError(f"Assembly summary file not found: {assembly}")

    # Deacon index: when provided for host depletion, must exist
    deacon_idx = args.get("deacon_index", "NA")
    if deacon_idx and str(deacon_idx).strip() not in ("", "NA"):
        try:
            validate_file_exists(deacon_idx, "Deacon index file")
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Deacon index file does not exist: {e}")


def get_samples_from_args(args: Dict[str, Any]) -> Dict[str, List[str]]:
    """Extract and validate samples from arguments.
    
    Args:
        args: Dictionary of pipeline arguments
        
    Returns:
        Dictionary mapping sample names to file paths
        
    Raises:
        ValidationError: If samples cannot be determined from arguments
    """
    sample_sheet = args.get("sample_sheet")
    samples = args.get("samples")
    data_type = args.get("data_type")
    
    if sample_sheet and os.path.isfile(sample_sheet):
        return validate_sample_sheet(sample_sheet, data_type)
    elif samples:
        return samples
    else:
        raise SampleConfigurationNotFoundError(
            "Either 'sample_sheet' or 'samples' must be provided. "
            f"Sample sheet path: {sample_sheet}"
        )

