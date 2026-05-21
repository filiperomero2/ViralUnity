#!/usr/bin/env python

"""
This scripts dynamically creates config files and runs
the viralunity consensus snakemake pipeline.
Filipe Moreira - 2024/09/21
"""

import logging
import os
import sys
from typing import Any, Dict

from viralunity import _orchestrator
from viralunity.constants import DataType, ResourceDefaults
from viralunity.validators import (
    CONSENSUS_PATH_ARG_KEYS,
    get_samples_from_args,
    resolve_path_args,
    validate_consensus_requirements,
    validate_illumina_requirements,
)

# Set up logging
logger = logging.getLogger(__name__)


def validate_args(args: Dict[str, Any]) -> Dict[str, list]:
    """Validate all pipeline arguments.

    Args:
        args: Dictionary of pipeline arguments

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        ValidationError: If validation fails
    """
    logger.info("Validating pipeline arguments")

    # Get and validate samples
    samples = get_samples_from_args(args)

    logger.info(f"Found {len(samples)} samples")

    # Validate consensus-specific requirements
    validate_consensus_requirements(args)

    # Handle primer scheme - set to "NA" if not provided
    if not args.get("primer_scheme"):
        logger.info("A primer scheme was not provided (untargeted sequencing).")
        args["primer_scheme"] = "NA"
    else:
        logger.info("A primer scheme was provided (Amplicon sequencing)...")

    # Validate data-type specific requirements
    validate_illumina_requirements(args)

    # ``validate_consensus_requirements`` parses ``--segmented-reference``
    # (``L=/path/L.fasta``) into a dict stored under ``reference``. Now
    # that the dict exists, resolve its values to absolute paths.
    if isinstance(args.get("reference"), dict):
        resolve_path_args(args, ("reference",))

    logger.info("All arguments validated successfully")
    return samples


def generate_config_file(samples: Dict[str, list], args: Dict[str, Any]) -> None:
    """Generate configuration file for the pipeline.

    Args:
        samples: Dictionary mapping sample names to file paths
        args: Dictionary of pipeline arguments

    Raises:
        ValidationError: If config generation fails
    """
    data_type = args["data_type"]
    generator = _orchestrator.start_config(args, samples)

    # Add consensus-specific settings
    generator.add_consensus_settings(
        reference=args["reference"],
        primer_scheme=args.get("primer_scheme", "NA"),
        minimum_coverage=args.get("minimum_coverage", 20),
        minimap2_consensus_align_flags=args.get(
            "minimap2_consensus_align_flags",
            "-a --sam-hit-only --secondary=no --score-N=0",
        ),
    )

    # Add workflow_path (consensus-specific)
    generator.add_workflow_path(sys.path[0])

    # Add Nanopore-specific settings if needed
    if data_type == DataType.NANOPORE:
        generator.add_consensus_nanopore_settings(
            minimum_read_length=args.get("minimum_read_length", 50),
            af_threshold=args.get("af_threshold", 0.51),
            chunk_size=args.get("chunk_size", 10000),
            clair3_model=args.get("clair3_model", "r1041_e82_400bps_sup_v500"),
            variant_quality=args.get("variant_quality", 20),
            variant_depth=args.get("variant_depth", 10),
            minimum_map_quality=args.get("minimum_map_quality", 30),
        )

    # Add Illumina-specific settings if needed
    if data_type == DataType.ILLUMINA:
        generator.add_illumina_settings(
            adapters=args["adapters"],
            minimum_read_length=args.get("minimum_read_length", 50),
            trim_head=args.get("trim_head", 0),
            trim_tail=args.get("trim_tail", 0),
            cut_front_mean_quality=args.get("cut_front_mean_quality", 10),
            cut_tail_mean_quality=args.get("cut_tail_mean_quality", 10),
            cut_right_window_size=args.get("cut_right_window_size", 4),
            cut_right_mean_quality=args.get("cut_right_mean_quality", 15),
            af_threshold=args.get("af_threshold", 0.51),
            af_isnv_threshold=args.get("af_isnv_threshold", 0),
            run_isnv=args.get("run_isnv", False),
        )

    # Add resource settings
    if data_type == DataType.ILLUMINA:
        generator.add_resource_settings(args, ResourceDefaults.CONSENSUS_ILLUMINA_RULES)
    else:
        generator.add_resource_settings(args, ResourceDefaults.CONSENSUS_NANOPORE_RULES)

    # Save config file
    generator.save()

    logger.info(f"Configuration file generated: {args['config_file']}")


def run_snakemake_workflow(args: Dict[str, Any]) -> bool:
    """Run the Snakemake workflow for the consensus pipeline."""
    thisdir = os.path.abspath(os.path.dirname(__file__))
    segmented_suffix = "_segmented" if isinstance(args.get("reference"), dict) else ""
    workflow_path = os.path.join(
        thisdir, "scripts", f"consensus_{args['data_type']}{segmented_suffix}.smk"
    )
    return _orchestrator.run_workflow(workflow_path, args)


def main(args: Dict[str, Any]) -> int:
    """Main entry point for the consensus pipeline.

    Returns:
        Exit code (0 for success, 1 for failure).
    """
    return _orchestrator.run_pipeline(
        args,
        resolve_paths=lambda a: resolve_path_args(a, CONSENSUS_PATH_ARG_KEYS),
        validate=validate_args,
        generate_config=generate_config_file,
        run_workflow_fn=run_snakemake_workflow,
    )
