#!/usr/bin/env python

"""
This scripts dynamically creates config files and runs
the viralunity metagenomics snakemake pipeline.
Filipe Moreira - 2024/09/21
"""

import os
import logging
from typing import Dict, Any
from pathlib import Path
from snakemake import snakemake

from viralunity.validators import (
    get_samples_from_args,
    validate_illumina_requirements,
    validate_nanopore_requirements,
    validate_metagenomics_requirements,
    resolve_path_args,
    META_PATH_ARG_KEYS,
)
from viralunity.config_generator import ConfigGenerator
from viralunity.exceptions import ValidationError
from viralunity.constants import DataType, ResourceDefaults

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
    data_type = args.get("data_type")
    if data_type == DataType.ILLUMINA:
        validate_illumina_requirements(args)
    elif data_type == DataType.NANOPORE:
        validate_nanopore_requirements(args)

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
    negative = args.get("negative_controls")
    if isinstance(negative, str):
        negative = [x.strip() for x in negative.split(",") if x.strip()]
    elif negative is None:
        negative = []

    # Normalize Krona database path: if it points to a dir containing 'taxonomy', append it.
    krona_db = args.get("krona_database", "NA")
    if krona_db and krona_db != "NA":
        kp = Path(krona_db)
        if (kp / "taxonomy").is_dir() and kp.name != "taxonomy":
            krona_db = str(kp / "taxonomy")

    generator.add_metagenomics_settings(
        kraken2_database=args["kraken2_database"],
        krona_database=krona_db,
        remove_human_reads=args.get("remove_human_reads", False),
        remove_unclassified_reads=args.get("remove_unclassified_reads", False),
        host_reference=args.get("host_reference", "NA"),
        deacon_index=args.get("deacon_index", "NA"),
        taxdump=args.get("taxdump", "NA"),
        run_denovo_assembly=args.get("run_denovo_assembly", False),
        run_kraken2_reads=args.get("run_kraken2_reads", True),
        run_kraken2_contigs=args.get("run_kraken2_contigs", True),
        run_diamond_reads=args.get("run_diamond_reads", False),
        run_diamond_contigs=args.get("run_diamond_contigs", False),
        taxids=args.get("taxids", "NA"),
        diamond_database=args.get("diamond_database", "NA"),
        diamond_sensitivity=args.get("diamond_sensitivity", "sensitive"),
        evalue=args.get("evalue", 0.001),
        bleed_fraction=args.get("bleed_fraction", 0.005),
        negative_controls=negative,
        negative_p_threshold=args.get("negative_p_threshold", 0.01),
        minimum_hit_group=args.get("minimum_hit_group", 4),
    )

    generator.add_reference_assembly_settings(
        run_reference_assembly=args.get("run_reference_assembly", False),
        method=args.get("method", "kraken2"),
        source=args.get("source", "reads"),
        reads_count=args.get("reads_count", 100),
        contigs_count=args.get("contigs_count", 1),
        families=args.get(
            "families",
            "Coronaviridae,Orthomyxoviridae,Flaviviridae,Herpesviridae,Papillomaviridae,Paramyxoviridae,Adenoviridae",
        ),
        reference_selection_strategy=args.get("reference_selection_strategy", "taxid"),
        blast_qcov=args.get("blast_qcov", 80),
        blast_pident=args.get("blast_pident", 80),
        viral_genomes=args.get("viral_genomes", "databases/virus_genomes/viral.genomes.fasta"),
        viral_taxids=args.get("viral_taxids", "databases/virus_genomes/genome2taxid.tsv"),
    )

    if data_type == DataType.ILLUMINA:
        generator.add_illumina_settings(
            adapters=args.get("adapters") or "NA",
            minimum_read_length=args.get("minimum_read_length", 50),
            trim_head=args.get("trim_head", 0),
            trim_tail=args.get("trim_tail", 0),
        )

    # Add resource settings
    shared_rules = ResourceDefaults.META_SHARED_RULES
    if data_type == DataType.ILLUMINA:
        generator.add_resource_settings(
            args, shared_rules + ResourceDefaults.META_ILLUMINA_RULES
        )
    else:
        generator.add_resource_settings(
            args, shared_rules + ResourceDefaults.META_NANOPORE_RULES
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
        thisdir, "scripts", f"metagenomics_{args['data_type']}.smk"
    )

    if not os.path.isfile(workflow_path):
        raise ValidationError(f"Workflow file not found: {workflow_path}")

    successful = snakemake(
        workflow_path,
        configfiles=[args["config_file"]],
        cores=args["threads_total"],
        use_conda=True,
        targets=["all"],
        forceall=True,
        lock=False,
        workdir=os.getcwd(),
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
        # Make every user-supplied path absolute relative to the shell's cwd
        # *before* validation or config generation looks at them. This means
        # the meaning of e.g. ``--host-reference databases/host/host.fasta``
        # no longer depends on where the user happens to place
        # ``--config-file``.
        resolve_path_args(args, META_PATH_ARG_KEYS)

        # Validate arguments
        samples = validate_args(args)

        if samples is None or len(samples) == 0:
            logger.warning("No samples were provided.")
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
