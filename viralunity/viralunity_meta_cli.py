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
        "--pipeline",
        type=str,
        help="Metagenomics pipeline version (default: v1). Use v2 for the extended nanopore workflow.",
        choices=["v1", "v2"],
        default="v1",
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
        help="Complete path for the kraken2 database directory (required for v1; optional for v2 if run_kraken2 is false)",
        required=False,
        default=None,
    )

    meta_parser.add_argument(
        "--krona-database",
        help="Complete path for the krona taxonomic database",
        required=True,
    )

    # v2 toggles (used by metagenomics_nanopore_v2 workflow)
    meta_parser.add_argument("--run-denovo-assembly", action="store_true", help="(v2) Run de novo assembly (MEGAHIT).")
    meta_parser.add_argument("--run-kraken2", action="store_true", help="(v2) Run Kraken2 on contigs.")
    meta_parser.add_argument("--run-diamond", action="store_true", help="(v2) Run DIAMOND blastx on contigs.")

    # DIAMOND parameters (v2)
    meta_parser.add_argument("--diamond-database", help="(v2) Path to DIAMOND protein FASTA database used to build .dmnd.")
    meta_parser.add_argument("--diamond-sensitivity", default="sensitive", help="(v2) DIAMOND sensitivity preset (e.g. sensitive, very-sensitive).")
    meta_parser.add_argument("--evalue", type=float, default=1e-10, help="(v2) DIAMOND e-value threshold.")

    # Host removal (v2)
    meta_parser.add_argument("--host-reference", default="NA", help="(v2) Host reference FASTA for dehosting; use NA to disable.")

    # Taxonomy resources used by v2 scripts
    meta_parser.add_argument("--taxdump", help="(v2) Path to NCBI taxdump directory.")
    meta_parser.add_argument("--assembly-summary", help="(v2) Path to viral_refseq assembly summary TSV.")
    meta_parser.add_argument("--taxid-to-family", help="(v2) Path to taxid_to_family CSV.")

    # Medaka
    #meta_parser.add_argument("--medaka-model", default="r941_min_high_g360", help="(v2) Medaka model name.")

    # Preferred (new) names
    meta_parser.add_argument(
        "--remove-human-sequences",
        action="store_true",
        help="Remove human-associated sequences (contigs or reads) from results",
    )

    meta_parser.add_argument(
        "--remove-unclassified-sequences",
        action="store_true",
        help="Remove unclassified sequences (taxid 0) from results",
    )

    # Backward-compatible aliases (hidden from help)
    meta_parser.add_argument(
        "--remove-human-reads",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    meta_parser.add_argument(
        "--remove-unclassified-reads",
        action="store_true",
        help=argparse.SUPPRESS,
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
