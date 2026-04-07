"""Click CLI for viralunity meta command."""

import click
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
        type=click.Choice(
            ["sensitive", "mid-sensitive", "more-sensitive", "ultra-sensitive"]
        ),
        default="sensitive",
        show_default=True,
        help="Diamond sensitivity mode.",
    ),
    click.option(
        "--evalue",
        default=0.001,
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
]


def _add_common_meta_options(func):
    """Decorator that stacks all shared meta options onto a click command."""
    for option in reversed(_COMMON_META_OPTIONS):
        func = option(func)
    return func


def _build_meta_args(data_type, **kwargs) -> dict:
    """Build args dict for meta_main, normalising negative_controls."""
    negative = kwargs.get("negative_controls", "")
    if isinstance(negative, str):
        kwargs["negative_controls"] = [
            x.strip() for x in negative.split(",") if x.strip()
        ]
    kwargs["data_type"] = data_type
    return kwargs


@click.group(name="meta")
def meta():
    """Metagenomics pipeline."""


@meta.command("illumina")
@_add_common_meta_options
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
def meta_illumina(
    sample_sheet,
    config_file,
    output,
    run_name,
    kraken2_database,
    krona_database,
    remove_human_reads,
    remove_unclassified_reads,
    host_reference,
    deacon_index,
    taxdump,
    run_denovo_assembly,
    run_kraken2_reads,
    run_kraken2_contigs,
    run_diamond_reads,
    run_diamond_contigs,
    taxids,
    diamond_database,
    diamond_sensitivity,
    evalue,
    bleed_fraction,
    negative_controls,
    negative_p_threshold,
    minimum_hit_group,
    threads,
    threads_total,
    create_config_only,
    adapters,
    minimum_read_length,
    trim_head,
    trim_tail,
    cut_front_mean_quality,
    cut_tail_mean_quality,
    cut_right_window_size,
    cut_right_mean_quality,
):
    """Run metagenomics pipeline for Illumina data."""
    args = _build_meta_args(
        data_type="illumina",
        sample_sheet=sample_sheet,
        config_file=config_file,
        output=output,
        run_name=run_name,
        kraken2_database=kraken2_database,
        krona_database=krona_database,
        remove_human_reads=remove_human_reads,
        remove_unclassified_reads=remove_unclassified_reads,
        host_reference=host_reference,
        deacon_index=deacon_index,
        taxdump=taxdump,
        run_denovo_assembly=run_denovo_assembly,
        run_kraken2_reads=run_kraken2_reads,
        run_kraken2_contigs=run_kraken2_contigs,
        run_diamond_reads=run_diamond_reads,
        run_diamond_contigs=run_diamond_contigs,
        taxids=taxids,
        diamond_database=diamond_database,
        diamond_sensitivity=diamond_sensitivity,
        evalue=evalue,
        bleed_fraction=bleed_fraction,
        negative_controls=negative_controls,
        negative_p_threshold=negative_p_threshold,
        minimum_hit_group=minimum_hit_group,
        threads=threads,
        threads_total=threads_total,
        create_config_only=create_config_only,
        adapters=adapters,
        minimum_read_length=minimum_read_length,
        trim_head=trim_head,
        trim_tail=trim_tail,
        cut_front_mean_quality=cut_front_mean_quality,
        cut_tail_mean_quality=cut_tail_mean_quality,
        cut_right_window_size=cut_right_window_size,
        cut_right_mean_quality=cut_right_mean_quality,
    )
    raise SystemExit(meta_main(args))


@meta.command("nanopore")
@_add_common_meta_options
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
def meta_nanopore(
    sample_sheet,
    config_file,
    output,
    run_name,
    kraken2_database,
    krona_database,
    remove_human_reads,
    remove_unclassified_reads,
    host_reference,
    deacon_index,
    taxdump,
    run_denovo_assembly,
    run_kraken2_reads,
    run_kraken2_contigs,
    run_diamond_reads,
    run_diamond_contigs,
    taxids,
    diamond_database,
    diamond_sensitivity,
    evalue,
    bleed_fraction,
    negative_controls,
    negative_p_threshold,
    minimum_hit_group,
    threads,
    threads_total,
    create_config_only,
    run_polish_racon,
    run_polish_medaka,
    medaka_model,
):
    """Run metagenomics pipeline for Nanopore data."""
    args = _build_meta_args(
        data_type="nanopore",
        sample_sheet=sample_sheet,
        config_file=config_file,
        output=output,
        run_name=run_name,
        kraken2_database=kraken2_database,
        krona_database=krona_database,
        remove_human_reads=remove_human_reads,
        remove_unclassified_reads=remove_unclassified_reads,
        host_reference=host_reference,
        deacon_index=deacon_index,
        taxdump=taxdump,
        run_denovo_assembly=run_denovo_assembly,
        run_kraken2_reads=run_kraken2_reads,
        run_kraken2_contigs=run_kraken2_contigs,
        run_diamond_reads=run_diamond_reads,
        run_diamond_contigs=run_diamond_contigs,
        taxids=taxids,
        diamond_database=diamond_database,
        diamond_sensitivity=diamond_sensitivity,
        evalue=evalue,
        bleed_fraction=bleed_fraction,
        negative_controls=negative_controls,
        negative_p_threshold=negative_p_threshold,
        minimum_hit_group=minimum_hit_group,
        threads=threads,
        threads_total=threads_total,
        create_config_only=create_config_only,
        run_polish_racon=run_polish_racon,
        run_polish_medaka=run_polish_medaka,
        medaka_model=medaka_model,
    )
    raise SystemExit(meta_main(args))
