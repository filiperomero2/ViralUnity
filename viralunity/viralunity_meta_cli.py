import argparse
from viralunity.viralunity_meta import main as meta_main


def _add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments common to all pipeline versions."""
    parser.add_argument(
        "--data-type",
        type=str,
        help="Sequencing data type (illumina or nanopore).",
        choices=["illumina", "nanopore"],
        nargs="?",
        required=True,
    )

    parser.add_argument(
        "--pipeline",
        type=str,
        help="Metagenomics pipeline version (default: v1). Use v2 for the extended nanopore workflow.",
        choices=["v1", "v2"],
        default="v1",
    )

    parser.add_argument(
        "--sample-sheet",
        help="Complete path for a csv file with samples data paths and metadata",
        required=True,
    )

    parser.add_argument(
        "--config-file",
        help="Complete path for input (viralunity config) file to be created.",
        required=True,
    )

    parser.add_argument(
        "--output",
        help="Complete path for output directory to be created by viralunity.",
        required=True,
    )

    parser.add_argument(
        "--run-name",
        type=str,
        help="Name for the sequencing run (optional).",
        nargs="?",
        const=1,
        default="undefined",
    )
    
    parser.add_argument(
        "--trim",
        type=int,
        help="Number of bases to trim from the 5' end of reads (Default = 0) [Illumina QC]",
        nargs="?",
        const=1,
        default=0,
    )

    parser.add_argument(
        "--create-config-only",
        help="Only create config file, not running the workflow (boolean)",
        action="store_true",
    )

    parser.add_argument(
        "--threads",
        type=int,
        help="Number of available threads for individual tasks (Default = 1)",
        nargs="?",
        const=1,
        default=1,
    )

    parser.add_argument(
        "--threads-total",
        type=int,
        help="Number of available threads for the entire workflow (Default = 1)",
        nargs="?",
        const=1,
        default=1,
    )

    parser.add_argument(
        "--krona-database",
        help="Complete path for the krona taxonomic database",
        required=True,
    )

    parser.add_argument(
        "--kraken2-database",
        help="Complete path for the kraken2 database directory (required for v1; optional for v2 if run_kraken2 is false)",
        required=False,
        default=None,
    )

    parser.add_argument(
        "--negative-controls",
        type=str,
        default="",
        help="Comma-separated sample names to treat as negative controls (e.g. NEG01,NEG02).",
    )

    parser.add_argument(
        "--bleed-fraction",
        type=float,
        default=0.005,
        help="Fraction of max RPM used for bleed-through filtering (default: 0.005).",
    )

    parser.add_argument(
        "--negative-p-threshold",
        type=float,
        default=0.01,
        help="P-value threshold for negative control background filter (default: 0.01).",
    )

    parser.add_argument(
        "--remove-human-sequences",
        action="store_true",
        help="Remove human-associated sequences (contigs or reads) from results",
    )

    parser.add_argument(
        "--remove-unclassified-sequences",
        action="store_true",
        help="Remove unclassified sequences (taxid 0) from results",
    )

    # Backward-compatible aliases (hidden from help)
    parser.add_argument(
        "--remove-human-reads",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--remove-unclassified-reads",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--adapters",
        help="Complete path for adapters sequences in fasta format [Illumina QC]",
    )

    parser.add_argument(
        "--minimum-read-length",
        type=int,
        help="Minimum read length threshold (Default = 50) [Illumina QC]",
        nargs="?",
        const=1,
        default=50,
    )


def _add_v1_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments specific to v1 pipeline."""
    # v1-specific arguments can be added here in the future
    pass


def _add_v2_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments specific to v2 pipeline."""
    # v2 workflow toggles
    parser.add_argument(
        "--run-denovo-assembly",
        action="store_true",
        help="(v2) Run de novo assembly (MEGAHIT)."
    )

    parser.add_argument(
        "--run-kraken2-contigs",
        action="store_true",
        help="(v2) Run Kraken2 classification on assembled contigs (requires --run-denovo-assembly)."
    )

    parser.add_argument(
        "--run-diamond-contigs",
        action="store_true",
        help="(v2) Run DIAMOND blastx on assembled contigs (requires --run-denovo-assembly)."
    )

    parser.add_argument(
        "--run-kraken2-reads",
        action="store_true",
        help="(v2) Run Kraken2 classification directly on reads (host-filtered fastq)."
    )

    parser.add_argument(
        "--run-diamond-reads",
        action="store_true",
        help="(v2) Run DIAMOND classification directly on reads (host-filtered fastq)."
    )

    # Hidden backward-compatibility flags
    parser.add_argument(
        "--run-kraken2",
        action="store_true",
        help=argparse.SUPPRESS
    )

    parser.add_argument(
        "--run-diamond",
        action="store_true",
        help=argparse.SUPPRESS
    )

    # DIAMOND parameters
    parser.add_argument(
        "--diamond-database",
        help="(v2) Path to DIAMOND protein FASTA database used to build .dmnd."
    )

    parser.add_argument(
        "--diamond-sensitivity",
        default="sensitive",
        help="(v2) DIAMOND sensitivity preset (e.g. sensitive, very-sensitive)."
    )

    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-10,
        help="(v2) DIAMOND e-value threshold."
    )

    # Host removal
    parser.add_argument(
        "--host-reference",
        default="NA",
        help="(v2) Host reference FASTA for dehosting; use NA to disable."
    )

    # Taxonomy resources
    parser.add_argument(
        "--taxdump",
        help="(v2) Path to NCBI taxdump directory."
    )

    parser.add_argument(
        "--assembly-summary",
        help="(v2) Path to viral_refseq assembly summary TSV."
    )

    parser.add_argument(
        "--taxid-to-family",
        help="(v2) Path to taxid_to_family CSV."
    )

    # Medaka (commented out for future use)
    # parser.add_argument(
    #     "--medaka-model",
    #     default="r941_min_high_g360",
    #     help="(v2) Medaka model name."
    # )


def fill_arg_parser_meta(subparsers: argparse._SubParsersAction):
    """Configure argument parser for the meta subcommand."""
    meta_parser = subparsers.add_parser(
        'meta',
        description="A script to generate config files and run the viralunity metagenomics pipeline"
    )
    meta_parser.set_defaults(func=meta_main)
    
    # Add common arguments (used by both v1 and v2)
    _add_common_arguments(meta_parser)
    
    # Add v1-specific arguments
    _add_v1_arguments(meta_parser)
    
    # Add v2-specific arguments
    _add_v2_arguments(meta_parser)
