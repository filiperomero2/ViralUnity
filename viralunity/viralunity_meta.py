#!/usr/bin/env python

"""
This scripts dynamically creates config files and runs
the viralunity metagenomics snakemake pipeline.
Filipe Moreira - 2024/09/21
"""

import os
import logging
from typing import Dict, Any
from snakemake import snakemake

from viralunity.validators import (
    get_samples_from_args,
    validate_illumina_requirements,
    validate_metagenomics_requirements
)
from viralunity.config_generator import ConfigGenerator
from viralunity.exceptions import ValidationError
from viralunity.constants import DataType

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
    
    # Validate metagenomics-specific requirements
    validate_metagenomics_requirements(args)
    
    # Validate data-type specific requirements
    validate_illumina_requirements(args)
    
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
    logger.info("Generating configuration file")
    
    data_type = args["data_type"]
    run_name = args["run_name"]
    
    generator = ConfigGenerator(args["config_file"])
    
    # Add samples
    generator.add_samples(samples, data_type)
    
    # Add common settings
    generator.add_output(args["output"], run_name)
    generator.add_threads(args["threads"])
    
    # Add metagenomics-specific settings
    generator.add_metagenomics_settings(
        kraken2_database=args["kraken2_database"],
        krona_database=args["krona_database"],
        remove_human_reads=args.get("remove_human_reads", False),
        remove_unclassified_reads=args.get("remove_unclassified_reads", False)
    )
    
    # Add Illumina-specific settings if needed
    if data_type == DataType.ILLUMINA:
        generator.add_illumina_settings(
            adapters=args["adapters"],
            minimum_read_length=args.get("minimum_read_length", 50),
            trim=args.get("trim", 0)
        )
    
    # Save config file
    generator.save()
    
    logger.info(f"Configuration file generated: {args['config_file']}")


def run_snakemake_workflow(args: Dict[str, Any]) -> bool:
    """Run the Snakemake workflow.
    
    Args:
        args: Dictionary of pipeline arguments
        
    Returns:
        True if workflow completed successfully, False otherwise
    """
    logger.info("Starting Snakemake workflow")
    
    thisdir = os.path.abspath(os.path.dirname(__file__))
    workflow_path = os.path.join(
        thisdir,
        'scripts',
        f"metagenomics_{args['data_type']}.smk"
    )
    
    if not os.path.isfile(workflow_path):
        raise ValidationError(f"Workflow file not found: {workflow_path}")
    
    successful = snakemake(
        workflow_path,
        configfiles=[args["config_file"]],
        cores=args["threads_total"],
        targets=["all"],
        forceall=True,
        lock=False,
        workdir=os.path.dirname(args["config_file"]),
    )
    
    if successful:
        logger.info("Snakemake workflow completed successfully")
    else:
        logger.error("Snakemake workflow failed")
    
    return successful


def main(args: Dict[str, Any]) -> int:
    """Main entry point for the metagenomics pipeline.
    
    Args:
        args: Dictionary of pipeline arguments
        
    Returns:
        Exit code (0 for success, 1 for failure)
    """
    try:
        # Validate arguments
        samples = validate_args(args)
        
        if(samples is None or len(samples) == 0):
            print("No samples were provided.")
            return 0
        
        # Generate config file
        generate_config_file(samples, args)
        
        # If only creating config, exit early
        if args.get("create_config_only", False):
            logger.info("Config file created. Exiting without running workflow.")
            return 0
        
        # Run workflow
        successful = run_snakemake_workflow(args)
        return 0 if successful else 1
        
    except ValidationError as e:
        logger.error(f"Validation error: {e}")
        return 1
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        return 1

