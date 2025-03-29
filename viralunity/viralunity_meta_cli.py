import argparse
from viralunity.viralunity_meta import main as meta_main

def fill_arg_parser_meta(subparsers: argparse._SubParsersAction):
    meta_parser = subparsers.add_parser('meta', description="A script to generate config files and run the viralunity metagenomics pipeline")
    meta_parser.set_defaults(func=meta_main)
    
    ## Add common arguments
    meta_parser.add_argument(
        "--data-type",
        type=str,
        help="Sequencing data type (illumina or nanopore).",
        choices=["illumina", "nanopore"],
        nargs="?",
        required=True,
    )

    meta_parser.add_argument(
        "--sample-sheet",
        help="Complete path for a csv file with samples data paths and metadata",
        required=True,
    )

    meta_parser.add_argument(
        "--config-file",
        help="Complete path for input (viralunity config) file to be created.",
        required=True,
    )

    meta_parser.add_argument(
        "--output",
        help="Complete path for output directory to be created by viralunity.",
        required=True,
    )

    meta_parser.add_argument(
        "--run-name",
        type=str,
        help="Name for the sequencing run (optional).",
        nargs="?",
        const=1,
        default="undefined",
    )
    
    meta_parser.add_argument(
        "--trim",
        type=int,
        help="Number of bases to trim from the 5' end of reads (Default = 0) [Illumina QC]",
        nargs="?",
        const=1,
        default=0,
    )

    meta_parser.add_argument(
        "--create-config-only",
        help="Only create config file, not running the workflow (boolean)",
        action="store_true",
    )

    meta_parser.add_argument(
        "--threads",
        type=int,
        help="Number of available threads for individual tasks (Default = 1)",
        nargs="?",
        const=1,
        default=1,
    )

    meta_parser.add_argument(
        "--threads-total",
        type=int,
        help="Number of available threads for the entire workflow (Default = 1)",
        nargs="?",
        const=1,
        default=1,
    )
    
    ## Add meta specific arguments
    meta_parser.add_argument(
        "--kraken2-database",
        help="Complete path for the kraken2 database directory",
        required=True,
    )

    meta_parser.add_argument(
        "--krona-database",
        help="Complete path for the krona taxonomic database",
        required=True,
    )

    meta_parser.add_argument(
        "--remove-human-reads",
        help="Remove human reads from krona plot (boolean)",
        action="store_true",
    )

    meta_parser.add_argument(
        "--remove-unclassified-reads",
        help="Remove unclassified reads from krona plot (boolean)",
        action="store_true",
    )

    meta_parser.add_argument(
        "--adapters",
        help="Complete path for adapters sequences in fasta format [Illumina QC]",
    )

    meta_parser.add_argument(
        "--minimum-read-length",
        type=int,
        help="Minimum read length threshold (Default = 50) [Illumina QC]",
        nargs="?",
        const=1,
        default=50,
    )
