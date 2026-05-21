"""Click CLI for viralunity meta command."""

from typing import Any

import click

from viralunity.constants import ResourceDefaults
from viralunity.viralunity_meta import main as meta_main

# Common options (applied to both illumina and nanopore subcommands)
_COMMON_META_OPTIONS = [
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
        "--kraken2-database",
        default="NA",
        show_default=True,
        help="Path to Kraken2 database.",
    ),
    click.option(
        "--krona-database",
        default="NA",
        show_default=True,
        help="Path to Krona taxonomic database.",
    ),
    click.option(
        "--remove-human-reads",
        is_flag=True,
        default=False,
        help="Remove human reads from Krona plot.",
    ),
    click.option(
        "--remove-unclassified-reads",
        is_flag=True,
        default=False,
        help="Remove unclassified reads from Krona plot.",
    ),
    click.option(
        "--host-reference",
        default="NA",
        show_default=True,
        help="Path to host genome FASTA for dehosting (minimap2). NA to skip.",
    ),
    click.option(
        "--deacon-index",
        default="NA",
        show_default=True,
        help="Path to Deacon minimizer index for host depletion (overrides --host-reference).",
    ),
    click.option(
        "--taxdump",
        default="NA",
        show_default=True,
        help="Path to NCBI taxdump directory (nodes.dmp, names.dmp).",
    ),
    click.option(
        "--run-denovo-assembly",
        is_flag=True,
        default=False,
        help="Run de novo assembly with MEGAHIT.",
    ),
    click.option(
        "--run-kraken2-reads/--no-kraken2-reads",
        default=True,
        show_default=True,
        help="Enable/disable Kraken2 classification of reads.",
    ),
    click.option(
        "--run-kraken2-contigs/--no-kraken2-contigs",
        default=True,
        show_default=True,
        help="Enable/disable Kraken2 classification of contigs.",
    ),
    click.option(
        "--run-diamond-reads/--no-diamond-reads",
        default=False,
        show_default=True,
        help="Enable/disable DIAMOND blastx on reads.",
    ),
    click.option(
        "--run-diamond-contigs/--no-diamond-contigs",
        default=False,
        show_default=True,
        help="Enable/disable DIAMOND on assembled contigs.",
    ),
    click.option(
        "--taxids",
        default="NA",
        show_default=True,
        help="protein2taxid.tsv mapping file produced by 'viralunity get-databases diamond'.",
    ),
    click.option(
        "--diamond-database",
        default="NA",
        show_default=True,
        help="Protein FASTA for Diamond database.",
    ),
    click.option(
        "--diamond-sensitivity",
        type=click.Choice(["sensitive", "mid-sensitive", "more-sensitive", "ultra-sensitive"]),
        default="sensitive",
        show_default=True,
        help="Diamond sensitivity mode.",
    ),
    click.option(
        "--evalue",
        default=0.0000000001,
        show_default=True,
        type=float,
        help="Diamond E-value threshold.",
    ),
    click.option(
        "--bleed-fraction",
        default=0.005,
        show_default=True,
        type=float,
        help="Max-RPM bleed filter fraction.",
    ),
    click.option(
        "--negative-controls",
        default="",
        show_default=True,
        help="Comma-separated sample IDs to use as negative controls.",
    ),
    click.option(
        "--negative-p-threshold",
        default=0.01,
        show_default=True,
        type=float,
        help="p-value threshold for negative-control filter.",
    ),
    click.option(
        "--minimum-hit-group",
        default=4,
        show_default=True,
        type=int,
        help="Kraken2 minimum-hit-group parameter.",
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
    click.option(
        "--run-reference-assembly/--no-run-reference-assembly",
        default=False,
        show_default=True,
        help="Enable/disable reference assembly from filtered taxonomic hits.",
    ),
    click.option(
        "--viral-genomes",
        default="NA",
        show_default=True,
        help="Path to the viral genomes FASTA file for reference assembly BLAST database/sequences.",
    ),
    click.option(
        "--viral-taxids",
        default="NA",
        show_default=True,
        help="Path to the genome2taxid mapping TSV file for reference assembly.",
    ),
    click.option(
        "--method",
        type=click.Choice(["kraken2", "diamond", "both"]),
        default=None,
        help="Method to use for identifying potential reference genomes.",
    ),
    click.option(
        "--source",
        type=click.Choice(["reads", "contigs", "both"]),
        default=None,
        help="Source of taxonomy data to use for identifying reference genomes.",
    ),
    click.option(
        "--reads-count",
        default=100,
        show_default=True,
        type=int,
        help="Minimum number of reads assigned to a viral family to trigger reference assembly.",
    ),
    click.option(
        "--contigs-count",
        default=1,
        show_default=True,
        type=int,
        help="Minimum number of contigs assigned to a viral family to trigger reference assembly.",
    ),
    click.option(
        "--families",
        default="Coronaviridae,Orthomyxoviridae,Flaviviridae,Herpesviridae,Papillomaviridae,Paramyxoviridae,Adenoviridae",
        show_default=True,
        help="Comma-separated list of viral families to perform reference assembly.",
    ),
    click.option(
        "--reference-selection-strategy",
        type=click.Choice(["taxid", "similarity"]),
        default="taxid",
        show_default=True,
        help="Strategy to associate a reference genome to a taxonomic assignment.",
    ),
    click.option(
        "--blast-qcov",
        default=80,
        show_default=True,
        type=int,
        help="Minimum query coverage for BLAST similarity reference selection.",
    ),
    click.option(
        "--blast-pident",
        default=80,
        show_default=True,
        type=int,
        help="Minimum percent identity for BLAST similarity reference selection.",
    ),
]


def _add_common_meta_options(func):
    """Decorator that stacks all shared meta options onto a click command."""
    for option in reversed(_COMMON_META_OPTIONS):
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


def _build_meta_args(data_type: str, **kwargs: Any) -> dict:
    """Build args dict for meta_main, normalising negative_controls."""
    negative = kwargs.get("negative_controls", "")
    if isinstance(negative, str):
        kwargs["negative_controls"] = [x.strip() for x in negative.split(",") if x.strip()]
    kwargs["data_type"] = data_type
    return kwargs


@click.group(name="meta")
def meta() -> None:
    """Viral metagenomics pipeline.

    \b
    Sub-commands branch on input data type:
    * illumina  paired-end short reads
    * nanopore  long reads

    Both flavours run optional Kraken2 / DIAMOND classification on reads
    and/or contigs, taxonomic summaries, Krona plots, and an optional
    reference-assembly checkpoint. Run ``viralunity meta <data_type> --help``
    for the full option set.
    """


@meta.command("illumina")
@_add_common_meta_options
@_add_resource_options(ResourceDefaults.META_SHARED_RULES + ResourceDefaults.META_ILLUMINA_RULES)
@click.option(
    "--adapters",
    default="NA",
    show_default=True,
    help="Path to adapter sequences FASTA (or NA for auto-detection) [fastp].",
)
@click.option(
    "--minimum-read-length",
    default=50,
    show_default=True,
    type=int,
    help="Minimum read length after trimming [fastp --length_required].",
)
@click.option(
    "--trim-head",
    default=None,
    type=int,
    help="Bases to trim from 5'.",
)
@click.option("--trim-tail", default=None, type=int, help="Bases to trim from 3'.")
@click.option(
    "--cut-front-mean-quality",
    default=20,
    show_default=True,
    type=int,
    help="cut_front mean quality threshold [fastp].",
)
@click.option(
    "--cut-tail-mean-quality",
    default=20,
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
    default=20,
    show_default=True,
    type=int,
    help="cut_right mean quality threshold [fastp].",
)
def meta_illumina(**kwargs: Any) -> None:
    """Run metagenomics pipeline for Illumina paired-end data.

    Performs adapter trimming (fastp), optional host depletion (Deacon),
    de novo assembly (MEGAHIT), and Kraken2 / DIAMOND classification on reads
    and/or contigs. Toggle classification stages with ``--run-kraken2-reads``,
    ``--run-kraken2-contigs``, ``--run-diamond-reads``,
    ``--run-diamond-contigs``.
    """
    args = _build_meta_args(data_type="illumina", **kwargs)
    raise SystemExit(meta_main(args))


@meta.command("nanopore")
@_add_common_meta_options
@_add_resource_options(ResourceDefaults.META_SHARED_RULES + ResourceDefaults.META_NANOPORE_RULES)
@click.option(
    "--run-polish-racon/--no-polish-racon",
    default=False,
    show_default=True,
    help="Run Racon polishing on MEGAHIT assembly.",
)
@click.option(
    "--run-polish-medaka/--no-polish-medaka",
    default=False,
    show_default=True,
    help="Run Medaka polishing on assembly.",
)
@click.option(
    "--medaka-model",
    default=None,
    help="Medaka model name (e.g. r941_min_high_g360). Uses Medaka default if omitted.",
)
@click.option(
    "--clair3-model",
    default=None,
    show_default=True,
    help=(
        "Clair3 model name used by reference-assembly consensus calling. "
        "Only used when --run-reference-assembly is set. "
        "Falls back to r1041_e82_400bps_sup_v500 if omitted."
    ),
)
def meta_nanopore(**kwargs: Any) -> None:
    """Run metagenomics pipeline for Nanopore long-read data.

    Performs optional host depletion (Deacon), de novo assembly (MEGAHIT),
    optional polishing (Racon / Medaka), and Kraken2 / DIAMOND classification
    on reads and/or contigs. Toggle polishing with ``--run-polish-racon`` /
    ``--run-polish-medaka``, select the Medaka model with ``--medaka-model``,
    and (when ``--run-reference-assembly`` is enabled) select the Clair3
    model with ``--clair3-model``.
    """
    args = _build_meta_args(data_type="nanopore", **kwargs)
    raise SystemExit(meta_main(args))
