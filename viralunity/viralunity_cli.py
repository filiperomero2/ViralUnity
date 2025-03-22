import argparse
import sys
from viralunity import __program__, __version__
from viralunity.viralunity_meta import main as meta_main
from viralunity.viralunity_consensus import main as consensus_main

def cli(args):
    parser = argparse.ArgumentParser(
        description="ViralUnity is a simple tool to perform analysis of viral high-throughput sequencing data.",
    )

    parser.add_argument('-v', '--version', action='version', version=f'{__program__} {__version__}')
    
    subparsers = parser.add_subparsers()
    parent = get_common_arg_parser()
    meta_parser = subparsers.add_parser('meta', description="A script to generate config files and run the viralunity metagenomics pipeline", parents=[parent])
    meta_parser.set_defaults(func=meta_main)
    consensus_parser = subparsers.add_parser('consensus',  description="A script to generate config files and run the viralunity consensus pipeline", parents=[parent])
    consensus_parser.set_defaults(func=consensus_main)
    
    fill_arg_parser_meta(meta_parser)
    fill_arg_parser_consensus(consensus_parser)
    
    arguments = parser.parse_args(args)
    try:
        arguments.func(vars(arguments))
    except Exception as e:
        print(f"Error: {e}")
        parser.print_help()
        raise e
    

def main():
    sys.exit(cli(sys.argv[1:]))
    
def get_common_arg_parser():
    parser = argparse.ArgumentParser(add_help=False
    )
    
    parser.add_argument(
        "--data-type",
        type=str,
        help="Sequencing data type (illumina or nanopore).",
        choices=["illumina", "nanopore"],
        nargs="?",
        required=True,
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
    
    return parser
    
def fill_arg_parser_consensus(parser):

    parser.add_argument(
        "--reference",
        help="Complete path for the reference genome in fasta format",
        required=True,
    )

    parser.add_argument(
        "--primer-scheme",
        help="Complete path for the primer scheme bed file (amplicon sequencing only).",
    )

    parser.add_argument(
        "--minimum-coverage",
        type=int,
        help="Minimum sequencing coverage for including base in consensus sequence (Default = 20)",
        nargs="?",
        const=1,
        default=20,
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
 
def fill_arg_parser_meta(parser):
    parser.add_argument(
        "--kraken2-database",
        help="Complete path for the kraken2 database directory",
        required=True,
    )

    parser.add_argument(
        "--krona-database",
        help="Complete path for the krona taxonomic database",
        required=True,
    )

    parser.add_argument(
        "--remove-human-reads",
        help="Remove human reads from krona plot (boolean)",
        action="store_true",
    )

    parser.add_argument(
        "--remove-unclassified-reads",
        help="Remove unclassified reads from krona plot (boolean)",
        action="store_true",
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

if __name__ == "__main__":
    main()