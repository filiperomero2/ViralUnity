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
        help="Bases to trim from 5' end of both reads (default: 0). Used as trim_head when --trim-head not set.",
        nargs="?",
        const=0,
        default=0,
    )

    meta_parser.add_argument(
        "--trim-head",
        type=int,
        default=None,
        help="Bases to trim from 5' end (overrides --trim for fastp).",
    )

    meta_parser.add_argument(
        "--trim-tail",
        type=int,
        default=None,
        help="Bases to trim from 3' end (default: 0).",
    )

    meta_parser.add_argument(
        "--cut-front-mean-quality",
        type=int,
        default=20,
        help="fastp cut_front mean quality threshold (default: 20).",
    )

    meta_parser.add_argument(
        "--cut-tail-mean-quality",
        type=int,
        default=20,
        help="fastp cut_tail mean quality threshold (default: 20).",
    )

    meta_parser.add_argument(
        "--cut-right-window-size",
        type=int,
        default=4,
        help="fastp cut_right window size (default: 4).",
    )

    meta_parser.add_argument(
        "--cut-right-mean-quality",
        type=int,
        default=20,
        help="fastp cut_right mean quality threshold (default: 20).",
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
    
    ## Add meta specific arguments (Kraken2 and/or Diamond; choose one or both)
    meta_parser.add_argument(
        "--kraken2-database",
        help="Path to Kraken2 database (required only if using Kraken2).",
        default="NA",
    )

    meta_parser.add_argument(
        "--krona-database",
        help="Path to Krona taxonomic database (required when running any classification).",
        default="NA",
    )

    meta_parser.add_argument(
        "--no-kraken2-reads",
        dest="run_kraken2_reads",
        action="store_false",
        help="Disable Kraken2 classification of reads (default: enabled).",
    )

    meta_parser.add_argument(
        "--run-kraken2-reads",
        dest="run_kraken2_reads",
        action="store_true",
        help="Enable Kraken2 classification of reads (default: enabled).",
    )

    meta_parser.add_argument(
        "--no-kraken2-contigs",
        dest="run_kraken2_contigs",
        action="store_false",
        help="Disable Kraken2 classification of contigs when assembly is run (default: enabled).",
    )

    meta_parser.add_argument(
        "--run-kraken2-contigs",
        dest="run_kraken2_contigs",
        action="store_true",
        help="Enable Kraken2 classification of contigs when assembly is run (default: enabled).",
    )

    meta_parser.add_argument(
        "--run-diamond-reads",
        dest="run_diamond_reads",
        action="store_true",
        default=False,
        help="Enable DIAMOND blastx on reads (requires --diamond-database and --assembly-summary).",
    )

    meta_parser.add_argument(
        "--no-diamond-reads",
        dest="run_diamond_reads",
        action="store_false",
        help="Disable DIAMOND on reads (default: disabled).",
    )

    meta_parser.add_argument(
        "--run-diamond-contigs",
        dest="run_diamond_contigs",
        action="store_true",
        default=False,
        help="Enable DIAMOND on assembled contigs when --run-denovo-assembly is set.",
    )

    meta_parser.add_argument(
        "--no-diamond-contigs",
        dest="run_diamond_contigs",
        action="store_false",
        help="Disable DIAMOND on contigs (default: disabled).",
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
        help="Path to adapter FASTA for fastp (optional: omit for auto-detection with --detect_adapter_for_pe).",
        default="NA",
    )

    meta_parser.add_argument(
        "--minimum-read-length",
        type=int,
        help="Minimum read length after trimming (default: 50) [fastp --length_required]",
        nargs="?",
        const=50,
        default=50,
    )

    # Dehosting, assembly, taxdump
    meta_parser.add_argument(
        "--host-reference",
        help="Path to host genome FASTA for dehosting (minimap2). Omit or set to NA to skip.",
        default="NA",
    )

    meta_parser.add_argument(
        "--deacon-index",
        help="Path to Deacon minimizer index for host depletion (alternative to --host-reference). If set, host depletion uses Deacon instead of minimap2.",
        default="NA",
    )

    meta_parser.add_argument(
        "--taxdump",
        help="Path to NCBI taxdump directory (nodes.dmp, names.dmp) for taxonomic summaries (required when running any classification).",
        default="NA",
    )

    meta_parser.add_argument(
        "--run-denovo-assembly",
        help="Run de novo assembly with MEGAHIT and classify contigs.",
        action="store_true",
    )

    # Nanopore polishing (ignored when --data-type illumina)
    meta_parser.add_argument(
        "--run-polish-racon",
        dest="run_polish_racon",
        action="store_true",
        default=False,
        help="[Nanopore] Run Racon polishing on MEGAHIT assembly (default: disabled).",
    )

    meta_parser.add_argument(
        "--no-polish-racon",
        dest="run_polish_racon",
        action="store_false",
        help="[Nanopore] Disable Racon polishing (default).",
    )

    meta_parser.add_argument(
        "--run-polish-medaka",
        dest="run_polish_medaka",
        action="store_true",
        default=False,
        help="[Nanopore] Run Medaka polishing on assembly (default: disabled).",
    )

    meta_parser.add_argument(
        "--no-polish-medaka",
        dest="run_polish_medaka",
        action="store_false",
        help="[Nanopore] Disable Medaka polishing (default).",
    )

    meta_parser.add_argument(
        "--medaka-model",
        type=str,
        default=None,
        help="[Nanopore] Medaka model name (e.g. r941_min_high_g360). Omit to use Medaka default.",
    )

    meta_parser.add_argument(
        "--assembly-summary",
        help="NCBI assembly_summary.txt for Diamond taxonomy (required if running Diamond).",
        default="NA",
    )

    meta_parser.add_argument(
        "--diamond-database",
        help="Protein FASTA for Diamond database (required if --run-diamond-reads or --run-diamond-contigs).",
        default="NA",
    )

    meta_parser.add_argument(
        "--diamond-sensitivity",
        choices=["sensitive", "mid-sensitive", "more-sensitive", "ultra-sensitive"],
        default="sensitive",
        help="Diamond sensitivity (default: sensitive).",
    )

    meta_parser.add_argument(
        "--evalue",
        type=float,
        default=0.001,
        help="Diamond E-value threshold (default: 0.001).",
    )

    meta_parser.add_argument(
        "--bleed-fraction",
        type=float,
        default=0.005,
        help="Max-RPM bleed filter fraction (default: 0.005).",
    )

    meta_parser.add_argument(
        "--negative-controls",
        help="Comma-separated sample IDs to use as negative controls for background filter.",
        default="",
    )

    meta_parser.add_argument(
        "--negative-p-threshold",
        type=float,
        default=0.01,
        help="p-value threshold for negative-control filter (default: 0.01).",
    )
