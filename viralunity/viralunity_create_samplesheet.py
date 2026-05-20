#!/usr/bin/env python

"""
Refactored version of viralunity_create_samplesheet.py demonstrating improved code design.

Create sample sheet from directory structure
Filipe Moreira - 2023/09/16
"""

import glob
import logging
import os
from typing import Dict, List

import click

from viralunity.constants import SampleSheetPattern, SampleSheetSeparator
from viralunity.exceptions import FileNotFoundError, ValidationError

logger = logging.getLogger(__name__)


@click.command("create-samplesheet")
@click.option(
    "--input",
    "input_dir",
    required=True,
    help="Directory containing sequencing data (level 0) or subdirectories with sequencing data (level 1).",
)
@click.option(
    "--output",
    required=True,
    help="Output CSV file path.",
)
@click.option(
    "--separator",
    type=click.Choice(
        [
            SampleSheetSeparator.UNDERSCORE,
            SampleSheetSeparator.HYPHEN,
            SampleSheetSeparator.DOT,
        ]
    ),
    default=SampleSheetSeparator.HYPHEN,
    show_default=True,
    help="Separator character used to extract the sample name from file/directory names.",
)
@click.option(
    "--pattern",
    type=click.Choice([SampleSheetPattern.R1, SampleSheetPattern.BARCODE]),
    default=SampleSheetPattern.R1,
    show_default=True,
    help="Pattern string for unpaired reads file names (used at level 0).",
)
@click.option(
    "--level",
    type=click.Choice(["0", "1"]),
    default="1",
    show_default=True,
    help="Directory level to search for sequencing data: 0 = files in --input, 1 = files in subdirectories.",
)
def create_samplesheet(input_dir, output, separator, pattern, level):
    """Generate a sample-sheet CSV file from a sequencing run directory."""
    args = {
        "input": input_dir,
        "output": output,
        "separator": separator,
        "pattern": pattern,
        "level": int(level),
    }
    try:
        validate_args(args)
        generate_sample_sheet(args)
        logger.info(f"Sample sheet file generated: {output}")
    except (ValidationError, FileNotFoundError) as e:
        raise click.ClickException(str(e))


def validate_args(args: Dict[str, any]) -> None:
    """Validate command line arguments.

    Args:
        args: Dictionary of arguments

    Raises:
        ValidationError: If validation fails
        FileNotFoundError: If input directory doesn't exist
    """
    input_dir = args["input"]

    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    logger.info(f"Input data directory: {input_dir}")

    level = args["level"]
    if level == 1:
        logger.info("Data will be searched for in subdirectories")
    else:
        logger.info("Data will be searched for only in this directory")

    output_file = args["output"]
    if os.path.isfile(output_file):
        raise ValidationError(f"Output file already exists: {output_file}")


def extract_sample_name(path: str, separator: str) -> str:
    """Extract sample name from file or directory path.

    Args:
        path: File or directory path
        separator: Separator character used to split the name

    Returns:
        Sample name (first part before separator)
    """
    basename = os.path.basename(path)
    return basename.split(separator)[0]


def find_files_in_directory(directory: str, pattern: str = None) -> List[str]:
    """Find files in a directory, optionally matching a pattern.

    Args:
        directory: Directory to search
        pattern: Optional pattern to match (e.g., '*R1*')

    Returns:
        List of file paths
    """
    if pattern:
        search_pattern = os.path.join(directory, f"*{pattern}*")
    else:
        search_pattern = os.path.join(directory, "*")

    files = glob.glob(search_pattern)
    return [f for f in files if os.path.isfile(f)]


def validate_sample_files(file_paths: List[str], sample_name: str) -> None:
    """Validate that sample files are valid.

    Args:
        file_paths: List of file paths
        sample_name: Name of the sample

    Raises:
        ValidationError: If file count is invalid
    """
    num_files = len(file_paths)
    if num_files not in [1, 2]:
        raise ValidationError(
            f"Invalid number of files for sample {sample_name}: "
            f"expected 1 or 2, found {num_files}"
        )


def find_samples_level_1(input_dir: str, separator: str) -> Dict[str, List[str]]:
    """Find samples in subdirectories (level 1).

    Args:
        input_dir: Base input directory
        separator: Separator character for sample names

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        ValidationError: If sample files are invalid
    """
    samples = {}
    directories = sorted(glob.glob(os.path.join(input_dir, "*")))

    for directory in directories:
        if not os.path.isdir(directory):
            continue

        sample_name = extract_sample_name(directory, separator)
        file_paths = sorted(find_files_in_directory(directory))

        validate_sample_files(file_paths, sample_name)
        samples[sample_name] = file_paths

        logger.debug(f"Found {len(file_paths)} file(s) for sample {sample_name}")

    return samples


def find_samples_level_0(input_dir: str, separator: str, pattern: str) -> Dict[str, List[str]]:
    """Find samples in the base directory (level 0).

    Args:
        input_dir: Base input directory
        separator: Separator character for sample names
        pattern: Pattern to match files (e.g., 'R1')

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        ValidationError: If sample files are invalid
    """
    samples = {}
    pattern_files = sorted(find_files_in_directory(input_dir, pattern))

    for pattern_file in pattern_files:
        sample_name = extract_sample_name(pattern_file, separator)

        # Find all files matching this sample name
        sample_files = sorted(find_files_in_directory(input_dir, sample_name))

        validate_sample_files(sample_files, sample_name)
        samples[sample_name] = sample_files

        logger.debug(f"Found {len(sample_files)} file(s) for sample {sample_name}")

    return samples


def generate_sample_sheet(args: Dict[str, any]) -> None:
    """Generate sample sheet CSV file.

    Args:
        args: Dictionary of arguments

    Raises:
        ValidationError: If sample sheet generation fails
    """
    input_dir = args["input"]
    level = args["level"]
    separator = args["separator"]
    pattern = args["pattern"]
    output_file = args["output"]

    logger.info("Searching for sample files")

    if level == 1:
        samples = find_samples_level_1(input_dir, separator)
    else:
        samples = find_samples_level_0(input_dir, separator, pattern)

    if not samples:
        raise ValidationError("No samples found in input directory")

    logger.info(f"Found {len(samples)} samples")

    # Write sample sheet
    try:
        with open(output_file, "w") as f:
            for sample_name, file_paths in sorted(samples.items()):
                if len(file_paths) == 2:
                    line = f"{sample_name},{file_paths[0]},{file_paths[1]}\n"
                else:
                    line = f"{sample_name},{file_paths[0]}\n"
                f.write(line)
    except IOError as e:
        raise ValidationError(f"Failed to write sample sheet: {e}")

    logger.info(f"Sample sheet generated: {output_file}")


if __name__ == "__main__":
    create_samplesheet()
