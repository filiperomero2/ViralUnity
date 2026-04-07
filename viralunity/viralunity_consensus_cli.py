"""Click CLI for viralunity consensus command."""

import click
from typing import Optional, Tuple
from viralunity.viralunity_consensus import main as consensus_main


def _parse_segmented_reference(ctx, param, value: Tuple[str, ...]):
    """Parse SEGMENT=PATH pairs into a dict, or return None if empty."""
    if not value:
        return None
    result = {}
    for item in value:
        if "=" not in item:
            raise click.BadParameter(
                f"Expected SEGMENT=PATH format, got: {item!r}",
                param=param,
            )
        segment, path = item.split("=", 1)
        result[segment.strip()] = path.strip()
    return result


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
        callback=_parse_segmented_reference,
        is_eager=False,
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


@click.group(name="consensus")
def consensus():
    """Consensus genome assembly pipeline."""


@consensus.command("illumina")
@_add_common_options
@click.option(
    "--adapters", default=None, help="Path to adapter sequences FASTA [fastp QC]."
)
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
    sample_sheet,
    config_file,
    output,
    run_name,
    reference,
    segmented_reference,
    primer_scheme,
    minimum_coverage,
    minimum_read_length,
    threads,
    threads_total,
    create_config_only,
    adapters,
    trim_head,
    trim_tail,
    cut_front_mean_quality,
    cut_tail_mean_quality,
    cut_right_window_size,
    cut_right_mean_quality,
    af_threshold,
    af_isnv_threshold,
    run_isnv,
):
    """Run consensus pipeline for Illumina data."""
    args = dict(
        data_type="illumina",
        sample_sheet=sample_sheet,
        config_file=config_file,
        output=output,
        run_name=run_name,
        reference=reference,
        segmented_reference=segmented_reference,
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
    raise SystemExit(consensus_main(args))


@consensus.command("nanopore")
@_add_common_options
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
    sample_sheet,
    config_file,
    output,
    run_name,
    reference,
    segmented_reference,
    primer_scheme,
    minimum_coverage,
    minimum_read_length,
    threads,
    threads_total,
    create_config_only,
    af_threshold,
    chunk_size,
    clair3_model,
    variant_quality,
    variant_depth,
    minimum_map_quality,
):
    """Run consensus pipeline for Nanopore data."""
    args = dict(
        data_type="nanopore",
        sample_sheet=sample_sheet,
        config_file=config_file,
        output=output,
        run_name=run_name,
        reference=reference,
        segmented_reference=segmented_reference,
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
    raise SystemExit(consensus_main(args))
