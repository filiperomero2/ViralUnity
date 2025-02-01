#!/usr/bin/env python

"""
This scripts dynamically creates config files and run 
the viralunity metagenomics snakemake pipeline.
Filipe Moreira - 2024/09/21
"""

import os
import sys
import argparse
from snakemake import snakemake
from viralunity.viralunity_consensus import validate_sample_sheet, define_job_id


def get_args(args):

    parser = argparse.ArgumentParser(
        description="A script to generate config files and run the viralunity metagenomics pipeline",
        usage="""viralunity_meta.py [args]""",
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

    return vars(parser.parse_args(args))


def validate_args(args):

    if os.path.isfile(args["sample_sheet"]):
        sample_sheet = args["sample_sheet"]
        samples = validate_sample_sheet(sample_sheet, args)
    else:
        print("Sample sheet file does not exist, please verify.")
        exit()

    if os.path.isfile(args["config_file"]):
        print("Config file already exists. Please specify a new one.")
        exit()

    if os.path.isdir(args["output"]):
        print("Output directory already exists. Please specify a new one.")
        exit()

    if os.path.isdir(args["kraken2_database"]):
        print("kraken2 database directory exists.")
    else:
        print("kraken2 database directory does not exist, please verify.")
        exit()

    if os.path.isdir(args["krona_database"]):
        print("krona taxonomic database directory exists.")
    else:
        print("krona taxonomic database does not exist, please verify.")
        exit()

    if args["data_type"] == "illumina":
        if os.path.isfile(args["adapters"]):
            print("Illumina adapter sequences file exists.")
        else:
            print("Illumina adapter sequences file not found, please verify.")
            exit()

    print("All arguments were succesfully verified.")

    return samples


def generate_config_file(samples, args):
    data_type = args["data_type"]
    with open(args["config_file"], "w") as f:

        samples_string = "samples:\n"
        f.write(samples_string)
        # TODO: refactor for clarity
        for sample in samples.keys():
            if data_type == "illumina":
                my_sample = (
                    f"    sample-{sample}: {samples[sample][0]} {samples[sample][1]}\n"
                )
            else:
                my_sample = f"    sample-{sample}: {samples[sample][0]}\n"
            f.write(my_sample)

        data_line = f"data: {data_type}\n"
        f.write(data_line)

        kraken2_database_line = f"kraken2_database: {args['kraken2_database']}\n"
        f.write(kraken2_database_line)

        krona_database_line = f"krona_database: {args['krona_database']}\n"
        f.write(krona_database_line)

        threads_line = f"threads: {str(args['threads'])}\n"
        f.write(threads_line)

        output_line = f"output: {args['output']}/{define_job_id(args)}/\n"
        f.write(output_line)

        # TODO: Analyze if it is possible to make an adapter pattern with illumina/nanopore
        if data_type == "illumina":
            adapters_line = f"adapters: {args['adapters']}\n"
            f.write(adapters_line)

            minimum_length_line = (
                f"minimum_length: {str(args['minimum_read_length'])}\n"
            )
            f.write(minimum_length_line)

            trim_line = f"trim: {str(args['trim'])}\n"
            f.write(trim_line)

        remove_human_reads_line = f"remove_human_reads: {args['remove_human_reads']}\n"
        f.write(remove_human_reads_line)

        remove_unclassified_reads_line = (
            f"remove_unclassified_reads: {args['remove_unclassified_reads']}\n"
        )
        f.write(remove_unclassified_reads_line)

        print("The config file was generated. -> ", args["config_file"])

    return


def main():
    args = get_args(sys.argv[1:])
    samples = validate_args(args)
    generate_config_file(samples, args)

    config = args["config_file"]
    cores = args["threads_total"]
    target_rule = "all"

    workflow_path = f"{sys.path[0]}/metagenomics_{args['data_type']}.smk"

    if args["create_config_only"]:
        print("Finished.")
        return 0

    successful = snakemake(
        workflow_path, configfiles=[config], cores=cores, targets=[target_rule]
    )

    return 0 if successful else 1


if __name__ == "__main__":
    sys.exit(main())
