"""Shared pipeline-orchestration helpers used by ``viralunity_consensus`` and
``viralunity_meta``.

Each pipeline still owns its own ``validate_args``, ``generate_config_file``,
``run_snakemake_workflow``, and ``main`` so that existing test patches at
those module-level names continue to work; the helpers in this module are
called from inside those functions to remove what would otherwise be
duplicated boilerplate.
"""

import logging
import os
from typing import Any, Callable, Dict, Optional

from snakemake import snakemake

from viralunity.config_generator import ConfigGenerator
from viralunity.exceptions import ValidationError

logger = logging.getLogger(__name__)


def start_config(args: Dict[str, Any], samples: Dict[str, list]) -> ConfigGenerator:
    """Instantiate a ``ConfigGenerator`` and add the keys both pipelines
    always set: samples, output, threads.

    The caller is responsible for adding pipeline-specific keys (consensus
    or metagenomics settings), per-rule resource settings, and calling
    ``generator.save()``.
    """
    logger.info("Generating configuration file")
    generator = ConfigGenerator(args["config_file"])
    generator.add_samples(samples, args["data_type"])
    generator.add_output(args["output"], args["run_name"])
    generator.add_threads(args["threads"])
    return generator


def run_workflow(workflow_path: str, args: Dict[str, Any]) -> bool:
    """Run a Snakemake workflow with the kwargs both pipelines share.

    Args:
        workflow_path: Absolute path to the ``.smk`` file to execute.
        args: Pipeline argument dict. Must contain ``config_file`` and
            ``threads_total``.

    Raises:
        ValidationError: If ``workflow_path`` does not exist.
    """
    logger.info("Starting Snakemake workflow")

    if not os.path.isfile(workflow_path):
        raise ValidationError(f"Workflow file not found: {workflow_path}")

    successful = snakemake(
        workflow_path,
        configfiles=[args["config_file"]],
        cores=args["threads_total"],
        use_conda=True,
        targets=["all"],
    )

    if successful:
        logger.info("Snakemake workflow completed successfully")
    else:
        logger.error("Snakemake workflow failed")

    return successful


def run_pipeline(
    args: Dict[str, Any],
    *,
    resolve_paths: Callable[[Dict[str, Any]], None],
    validate: Callable[[Dict[str, Any]], Optional[Dict[str, list]]],
    generate_config: Callable[[Dict[str, list], Dict[str, Any]], None],
    run_workflow_fn: Callable[[Dict[str, Any]], bool],
    skip_when_no_samples: bool = False,
) -> int:
    """The try/except ``main`` skeleton shared by both pipelines.

    Args:
        args: Pipeline argument dict.
        resolve_paths: Called first to convert relative path args to
            absolute paths.
        validate: Called to validate args and return the samples dict.
        generate_config: Called to write the Snakemake config file.
        run_workflow_fn: Called to execute the Snakemake workflow. The
            pipeline-specific module passes its own
            ``run_snakemake_workflow`` so test patches at that name
            keep working.
        skip_when_no_samples: If ``True``, return 0 (success) when the
            validated samples dict is empty. ``viralunity_meta`` opts in;
            ``viralunity_consensus`` does not.

    Returns:
        Exit code (0 for success, 1 for failure).
    """
    try:
        resolve_paths(args)
        samples = validate(args)

        if skip_when_no_samples and (samples is None or len(samples) == 0):
            logger.warning("No samples were provided.")
            return 0

        generate_config(samples, args)

        if args.get("create_config_only", False):
            logger.info("Config file created. Exiting without running workflow.")
            return 0

        successful = run_workflow_fn(args)
        return 0 if successful else 1

    except ValidationError as e:
        logger.error(f"Validation error: {e}")
        return 1
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        return 1
