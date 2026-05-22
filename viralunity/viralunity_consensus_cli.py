"""Click CLI for viralunity consensus command."""

from typing import Any, Optional, Tuple

import click

from viralunity.constants import ResourceDefaults
from viralunity.viralunity_consensus import main as consensus_main


def _parse_segmented_reference(segmented_reference: Tuple[str, ...]) -> Optional[dict]:
    """Parse SEGMENT=PATH strings into a dictionary.

    Args:
        segmented_reference: Tuple of strings in SEGMENT=PATH format

    Returns:
        Dictionary mapping segment names to paths, or None if empty
    """
    if not segmented_reference:
        return None

    parsed = {}
    for entry in segmented_reference:
        if "=" not in entry:
            raise click.BadParameter(
                f"Invalid segmented reference format: '{entry}'. "
                "Expected SEGMENT=PATH (e.g. S=/path/to/S.fasta)"
            )
        name, path = entry.split("=", 1)
        name = name.strip()
        path = path.strip()
        if not name or not path:
            raise click.BadParameter(
                f"Invalid segmented reference format: '{entry}'. "
                "Both segment name and path are required."
            )
        parsed[name] = path
    return parsed


# Common options (applied to both illumina and nanopore subcommands)
_COMMON_OPTIONS = [
    click.option(
        "--sample-sheet",
        required=True,
        help="Path to CSV file with sample data paths and metadata.",
    ),
    click.option(
        "--config-file",
        required=True,
        help="Path for the config YAML file to be created.",
    ),
    click.option("--output", required=True, help="Path for the output directory."),
    click.option(
        "--run-name",
        default="undefined",
        show_default=True,
        help="Name for the sequencing run.",
    ),
    click.option(
        "--reference",
        default=None,
        help="Path to reference genome FASTA (mutually exclusive with --segmented-reference).",
    ),
    click.option(
        "--segmented-reference",
        multiple=True,
        metavar="SEGMENT=PATH",
        help="Reference for segmented viruses: SEGMENT=PATH per segment "
        "(e.g. S=/path/S.fasta). Mutually exclusive with --reference.",
    ),
    click.option(
        "--primer-scheme",
        default=None,
        help="Path to primer scheme BED file (amplicon sequencing only).",
    ),
    click.option(
        "--minimum-coverage",
        default=20,
        show_default=True,
        type=int,
        help="Minimum coverage depth for consensus base inclusion.",
    ),
    click.option(
        "--minimum-read-length",
        default=50,
        show_default=True,
        type=int,
        help="Minimum read length threshold.",
    ),
    click.option(
        "--threads",
        default=1,
        show_default=True,
        type=int,
        help="Threads for individual tasks.",
    ),
    click.option(
        "--threads-total",
        default=1,
        show_default=True,
        type=int,
        help="Total threads for the entire workflow.",
    ),
    click.option(
        "--create-config-only",
        is_flag=True,
        default=False,
        help="Only create the config file; do not run the workflow.",
    ),
]


def _add_common_options(func):
    """Decorator that stacks all shared options onto a click command."""
    for option in reversed(_COMMON_OPTIONS):
        func = option(func)
    return func


def _generate_resource_options(rules: list) -> list:
    """Generate click options for CPUs and RAM for a list of rules."""
    options = []
    for rule in rules:
        cmd_rule = rule.replace("_", "-")
        options.append(
            click.option(
                f"--{cmd_rule}-cpus",
                default=ResourceDefaults.DEFAULT_CPUS,
                show_default=True,
                type=int,
                help=f"Threads for {rule} rule.",
            )
        )
        options.append(
            click.option(
                f"--{cmd_rule}-ram",
                default=ResourceDefaults.DEFAULT_RAM,
                show_default=True,
                type=int,
                help=f"RAM (GB) for {rule} rule.",
            )
        )
    return options


def _add_resource_options(rules: list):
    """Decorator to add resource options to a click command."""

    def decorator(func):
        for option in reversed(_generate_resource_options(rules)):
            func = option(func)
        return func

    return decorator


@click.group(name="consensus")
def consensus() -> None:
    """Reference-guided consensus genome assembly pipeline.

    \b
    Sub-commands branch on input data type:
    * illumina  paired-end short reads (fastp + minimap2 + LoFreq)
    * nanopore  long reads (minimap2 + Clair3)

    Run ``viralunity consensus <data_type> --help`` for the full option set.
    """


@consensus.command("illumina")
@_add_common_options
@_add_resource_options(ResourceDefaults.CONSENSUS_ILLUMINA_RULES)
@click.option("--adapters", default=None, help="Path to adapter sequences FASTA [fastp QC].")
@click.option(
    "--trim-head",
    default=0,
    show_default=True,
    type=int,
    help="Bases to trim from 5' end [fastp QC].",
)
@click.option(
    "--trim-tail",
    default=0,
    show_default=True,
    type=int,
    help="Bases to trim from 3' end [fastp QC].",
)
@click.option(
    "--cut-front-mean-quality",
    default=10,
    show_default=True,
    type=int,
    help="cut_front mean quality threshold [fastp].",
)
@click.option(
    "--cut-tail-mean-quality",
    default=10,
    show_default=True,
    type=int,
    help="cut_tail mean quality threshold [fastp].",
)
@click.option(
    "--cut-right-window-size",
    default=4,
    show_default=True,
    type=int,
    help="cut_right window size [fastp].",
)
@click.option(
    "--cut-right-mean-quality",
    default=15,
    show_default=True,
    type=int,
    help="cut_right mean quality threshold [fastp].",
)
@click.option(
    "--af-threshold",
    default=0.51,
    show_default=True,
    type=float,
    help="Minimum allele frequency to call a variant into consensus.",
)
@click.option(
    "--af-isnv-threshold",
    default=0.0,
    show_default=True,
    type=float,
    help="Minimum allele frequency to call a variant into iSNV analysis.",
)
@click.option(
    "--run-isnv",
    is_flag=True,
    default=False,
    help="Run intra-host SNV analysis with LoFreq.",
)
def consensus_illumina(
    sample_sheet: str,
    config_file: str,
    output: str,
    run_name: str,
    reference: Optional[str],
    segmented_reference: Tuple[str, ...],
    primer_scheme: Optional[str],
    minimum_coverage: int,
    minimum_read_length: int,
    threads: int,
    threads_total: int,
    create_config_only: bool,
    adapters: Optional[str],
    trim_head: int,
    trim_tail: int,
    cut_front_mean_quality: int,
    cut_tail_mean_quality: int,
    cut_right_window_size: int,
    cut_right_mean_quality: int,
    af_threshold: float,
    af_isnv_threshold: float,
    run_isnv: bool,
    **kwargs: Any,
) -> None:
    """Run consensus pipeline for Illumina paired-end data.

    Performs adapter trimming (fastp), reference alignment (minimap2),
    optional primer trimming (ivar), variant calling (LoFreq), and consensus
    generation. Enable intra-host SNV analysis with ``--run-isnv``.

    For segmented references (e.g. influenza), pass each segment as
    ``--segmented-reference SEGMENT=PATH``; otherwise use ``--reference``.
    """
    args = dict(
        data_type="illumina",
        sample_sheet=sample_sheet,
        config_file=config_file,
        output=output,
        run_name=run_name,
        reference=reference,
        segmented_reference=_parse_segmented_reference(segmented_reference),
        primer_scheme=primer_scheme,
        minimum_coverage=minimum_coverage,
        minimum_read_length=minimum_read_length,
        threads=threads,
        threads_total=threads_total,
        create_config_only=create_config_only,
        adapters=adapters,
        trim_head=trim_head,
        trim_tail=trim_tail,
        cut_front_mean_quality=cut_front_mean_quality,
        cut_tail_mean_quality=cut_tail_mean_quality,
        cut_right_window_size=cut_right_window_size,
        cut_right_mean_quality=cut_right_mean_quality,
        af_threshold=af_threshold,
        af_isnv_threshold=af_isnv_threshold,
        run_isnv=run_isnv,
    )
    args.update(kwargs)
    raise SystemExit(consensus_main(args))


@consensus.command("nanopore")
@_add_common_options
@_add_resource_options(ResourceDefaults.CONSENSUS_NANOPORE_RULES)
@click.option(
    "--af-threshold",
    default=0.51,
    show_default=True,
    type=float,
    help="Minimum allele frequency to call a variant into consensus.",
)
@click.option(
    "--chunk-size",
    default=10000,
    show_default=True,
    type=int,
    help="Chunk size for clair3 processing.",
)
@click.option(
    "--clair3-model",
    default="r1041_e82_400bps_sup_v500",
    show_default=True,
    help="Clair3 model for variant calling.",
)
@click.option(
    "--variant-quality",
    default=20,
    show_default=True,
    type=int,
    help="Minimum variant quality (clair3).",
)
@click.option(
    "--variant-depth",
    default=10,
    show_default=True,
    type=int,
    help="Minimum alt allele depth to call variant (clair3).",
)
@click.option(
    "--minimum-map-quality",
    default=30,
    show_default=True,
    type=int,
    help="Minimum mapping quality (clair3).",
)
def consensus_nanopore(
    sample_sheet: str,
    config_file: str,
    output: str,
    run_name: str,
    reference: Optional[str],
    segmented_reference: Tuple[str, ...],
    primer_scheme: Optional[str],
    minimum_coverage: int,
    minimum_read_length: int,
    threads: int,
    threads_total: int,
    create_config_only: bool,
    af_threshold: float,
    chunk_size: int,
    clair3_model: str,
    variant_quality: int,
    variant_depth: int,
    minimum_map_quality: int,
    **kwargs: Any,
) -> None:
    """Run consensus pipeline for Nanopore long-read data.

    Performs reference alignment (minimap2), variant calling (Clair3 with a
    user-selectable model via ``--clair3-model``), and consensus generation.
    Reads shorter than ``--minimum-read-length`` are filtered upstream.

    For segmented references (e.g. influenza), pass each segment as
    ``--segmented-reference SEGMENT=PATH``; otherwise use ``--reference``.
    """
    args = dict(
        data_type="nanopore",
        sample_sheet=sample_sheet,
        config_file=config_file,
        output=output,
        run_name=run_name,
        reference=reference,
        segmented_reference=_parse_segmented_reference(segmented_reference),
        primer_scheme=primer_scheme,
        minimum_coverage=minimum_coverage,
        minimum_read_length=minimum_read_length,
        threads=threads,
        threads_total=threads_total,
        create_config_only=create_config_only,
        af_threshold=af_threshold,
        chunk_size=chunk_size,
        clair3_model=clair3_model,
        variant_quality=variant_quality,
        variant_depth=variant_depth,
        minimum_map_quality=minimum_map_quality,
    )
    args.update(kwargs)
    raise SystemExit(consensus_main(args))
