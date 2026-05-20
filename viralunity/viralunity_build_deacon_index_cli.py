"""Click CLI for viralunity build-deacon-index command."""

import logging
from pathlib import Path

import click

from viralunity._subprocess import run_command

logger = logging.getLogger(__name__)


@click.command("build-deacon-index")
@click.option(
    "--path",
    default="databases",
    show_default=True,
    help="Parent directory where the deacon_indexes/ subdirectory will be created.",
)
@click.option(
    "--input",
    "input_fasta",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Input Fasta file to index.",
)
@click.option(
    "--threads",
    default=8,
    show_default=True,
    help="Number of threads to use.",
)
def build_deacon_index(path, input_fasta, threads):
    """Build a Deacon minimizer index from a FASTA file.

    Creates {path}/deacon_indexes/ and runs 'deacon index build' on the input file.
    """
    db_dir = Path(path) / "deacon_indexes"
    db_dir.mkdir(parents=True, exist_ok=True)

    input_path = Path(input_fasta)
    output_file = db_dir / f"{input_path.stem}.idx"

    cmd = [
        "deacon",
        "index",
        "build",
        str(input_path),
        "-o",
        str(output_file),
        "-t",
        str(threads),
    ]

    click.echo(f"Building Deacon index for '{input_path.name}'...")
    run_command(cmd)

    click.echo("\nDeacon index built successfully.")
    click.echo(f"  Index: {output_file}")
    click.echo(f"Use --deacon-index {output_file} in your viralunity meta commands.")
