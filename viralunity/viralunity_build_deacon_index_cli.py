"""Click CLI for viralunity build-deacon-index command."""

import logging
import shutil
import subprocess
from pathlib import Path

import click

logger = logging.getLogger(__name__)


def _run(cmd: list, cwd: str | None = None) -> None:
    """Run a shell command, streaming output and raising on failure."""
    click.echo(f"$ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        raise click.ClickException(
            f"Command failed with exit code {result.returncode}: {' '.join(str(c) for c in cmd)}"
        )


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
    _run(cmd)

    click.echo(f"\nDeacon index built successfully.")
    click.echo(f"  Index: {output_file}")
    click.echo(f"Use --deacon-index {output_file} in your viralunity meta commands.")
