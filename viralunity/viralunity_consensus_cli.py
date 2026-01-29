import argparse
from viralunity.viralunity_consensus import main as consensus_main

def fill_arg_parser_consensus(subparsers: argparse._SubParsersAction):
    consensus_parser = subparsers.add_parser('consensus',  description="A script to generate config files and run the viralunity consensus pipeline")
    consensus_parser.set_defaults(func=consensus_main)
    
    
    ## Add common arguments
    consensus_parser.add_argument(
        "--data-type",
        type=str,
        help="Sequencing data type (illumina or nanopore).",
        choices=["illumina", "nanopore"],
        nargs="?",
        required=True,
    )

    consensus_parser.add_argument(
        "--sample-sheet",
        help="Complete path for a csv file with samples data paths and metadata",
        required=True,
    )

    consensus_parser.add_argument(
        "--config-file",
        help="Complete path for input (viralunity config) file to be created.",
        required=True,
    )

    consensus_parser.add_argument(
        "--output",
        help="Complete path for output directory to be created by viralunity.",
        required=True,
    )

    consensus_parser.add_argument(
        "--run-name",
        type=str,
        help="Name for the sequencing run (optional).",
        nargs="?",
        const=1,
        default="undefined",
    )
    
    consensus_parser.add_argument(
        "--trim-head",
        type=int,
        help="Number of bases to trim from the 5' end of reads (Default = 0) [Illumina QC]",
        nargs="?",
        const=1,
        default=0,
    )

    consensus_parser.add_argument(
        "--trim-tail",
        type=int,
        help="Number of bases to trim from the 3' end of reads (Default = 0) [Illumina QC]",
        nargs="?",
        const=1,
        default=0,
    )

    consensus_parser.add_argument(
        "--cut-front-mean-quality",
        type=int,
        help="Mean quality requirement option for cut_front (Default = 10) [Illumina QC]",
        nargs="?",
        const=1,
        default=10,
    )

    consensus_parser.add_argument(
        "--cut-tail-mean-quality",
        type=int,
        help="Mean quality requirement option for cut_tail (Default = 10) [Illumina QC]",
        nargs="?",
        const=1,
        default=10,
    )

    consensus_parser.add_argument(
        "--cut-right-window-size",
        type=int,
        help="Window size for cut_right (Default = 4) [Illumina QC]",
        nargs="?",
        const=1,
        default=4,
    )

    consensus_parser.add_argument(
        "--cut-right-mean-quality",
        type=int,
        help="Mean quality requirement option for cut_right (Default = 15) [Illumina QC]",
        nargs="?",
        const=1,
        default=15,
    )

    consensus_parser.add_argument(
        "--af-threshold",
        type=float,
        help="Allele frequency threshold to call a variant into consensus (Default = 0.51)",
        nargs="?",
        const=1,
        default=0.51,
    )

    consensus_parser.add_argument(
        "--create-config-only",
        help="Only create config file, not running the workflow (boolean)",
        action="store_true",
    )

    consensus_parser.add_argument(
        "--threads",
        type=int,
        help="Number of available threads for individual tasks (Default = 1)",
        nargs="?",
        const=1,
        default=1,
    )

    consensus_parser.add_argument(
        "--threads-total",
        type=int,
        help="Number of available threads for the entire workflow (Default = 1)",
        nargs="?",
        const=1,
        default=1,
    )
    
    ## Add consensus specific arguments
    consensus_parser.add_argument(
        "--reference",
        help="Complete path for the reference genome in fasta format",
        required=True,
    )

    consensus_parser.add_argument(
        "--primer-scheme",
        help="Complete path for the primer scheme bed file (amplicon sequencing only).",
    )

    consensus_parser.add_argument(
        "--minimum-coverage",
        type=int,
        help="Minimum sequencing coverage for including base in consensus sequence (Default = 20)",
        nargs="?",
        const=1,
        default=20,
    )

    consensus_parser.add_argument(
        "--adapters",
        help="Complete path for adapters sequences in fasta format [Illumina QC]",
    )

    consensus_parser.add_argument(
        "--minimum-read-length",
        type=int,
        help="Minimum read length threshold (Default = 50)",
        nargs="?",
        const=1,
        default=50,
    )
 
