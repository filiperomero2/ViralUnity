"""Validation functions for ViralUnity pipeline arguments and data."""

import os
from typing import Dict, Any, List
import pandas as pd

from viralunity.viralunity.exceptions import (
    Kraken2DatabaseNotFoundError,
    KronaDatabaseNotFoundError,
    ValidationError,
    FileNotFoundError,
    SampleSheetError,
    SampleConfigurationNotFoundError,
    AdaptersNotFoundError,
    ReferenceNotFoundError,
    PrimerSchemeNotFoundError
)
from viralunity.viralunity.constants import DataType


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
    """Validate Illumina-specific requirements.
    
    Args:
        args: Dictionary of pipeline arguments
        
    Raises:
        ValidationError: If Illumina requirements are not met
    """
    if args.get("data_type") != DataType.ILLUMINA:
        return
    
    adapters = args.get("adapters")
    if not adapters:
        raise AdaptersNotFoundError("Illumina adapter sequences file is required")
    
    try:
        validate_file_exists(adapters, "Illumina adapter sequences file")
    except FileNotFoundError as e:
        raise AdaptersNotFoundError(f"Illumina adapter sequences file does not exist: {e}")


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
    
    Args:
        args: Dictionary of pipeline arguments
        
    Raises:
        ValidationError: If metagenomics requirements are not met
    """
    kraken2_db = args.get("kraken2_database")
    if not kraken2_db:
        raise Kraken2DatabaseNotFoundError("Kraken2 database directory is required")
    
    try:
        validate_directory_exists(kraken2_db, "Kraken2 database directory")
    except FileNotFoundError as e:
        raise Kraken2DatabaseNotFoundError(f"Kraken2 database directory does not exist: {e}")
    
    krona_db = args.get("krona_database")
    if not krona_db:
        raise KronaDatabaseNotFoundError("Krona database directory is required")
    
    try:
        validate_directory_exists(krona_db, "Krona database directory")
    except FileNotFoundError as e:
        raise KronaDatabaseNotFoundError(f"Krona database directory does not exist: {e}")


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

