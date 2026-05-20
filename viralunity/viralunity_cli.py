"""Top-level click CLI for ViralUnity."""

import click

from viralunity import __program__, __version__
from viralunity.viralunity_build_deacon_index_cli import build_deacon_index
from viralunity.viralunity_consensus_cli import consensus
from viralunity.viralunity_create_samplesheet import create_samplesheet
from viralunity.viralunity_get_databases_cli import get_databases
from viralunity.viralunity_meta_cli import meta


@click.group()
@click.version_option(version=__version__, prog_name=__program__)
def cli():
    """ViralUnity is a simple tool to perform analysis of viral high-throughput sequencing data."""


cli.add_command(consensus)
cli.add_command(meta)
cli.add_command(create_samplesheet)
cli.add_command(get_databases)
cli.add_command(build_deacon_index)


def main():
    cli()


if __name__ == "__main__":
    main()
