#!/usr/bin/env python

"""
This scripts dynamically creates config files and run 
the viralunity consensus snakemake pipeline.
Filipe Moreira - 2024/09/21
"""

import os
import sys
import pandas as pd
import datetime
from snakemake import snakemake

def validate_args(args):

    if os.path.isfile(args["sample_sheet"]):
        sample_sheet = args["sample_sheet"]
        samples = validate_sample_sheet(sample_sheet, args)
    else:
        raise Exception(f"Sample sheet file does not exist: {args['sample_sheet']}")

    if os.path.isfile(args["config_file"]):
        raise Exception("Config file already exists.")

    if os.path.isdir(args["output"]):
        raise Exception(f"Output directory '{args['output']}' already exists.")

    if not os.path.isfile(args["reference"]):
        raise Exception("Reference sequence file does not exist.")

    if args["primer_scheme"]:
        print("A primer scheme was provided (Amplicon sequencing)...")
        if not os.path.isfile(args["primer_scheme"]):
            raise Exception("Primer scheme file not found.")
    else:
        print("A primer scheme was not provided (untargeted sequencing).")
        args["primer_scheme"] = "NA"

    if args["data_type"] == "illumina":
        if args["adapters"]:
            print("Illumina adapter sequences were provided...")
            if not os.path.isfile(args["adapters"]):
                print("Adapter sequences file not found, please verify.")
                raise Exception("Illumina adapter sequences file not found.")
        else:
            raise Exception("Illumina adapter sequences file not found.")

    print("All arguments were succesfully verified.")

    return samples


def validate_sample_sheet(sample_sheet, args):
    df = pd.read_csv(sample_sheet, header=None)
    samples = {}
    for ind in df.index:
        if args["data_type"] == "illumina":
            # This assumes all illumina data are paired-end reads
            sample_name = df[0][ind]
            sample_R1 = df[1][ind]
            sample_R2 = df[2][ind]
            if os.path.isfile(sample_R1) and os.path.isfile(sample_R2):
                samples[sample_name] = [sample_R1, sample_R2]
                print("Files for sample", sample_name, "identified.")
            else:
                print("Problem in reads files for sample", sample_name)
                print("Please specify correct file paths.")
                raise Exception(f"Problem in reads files for sample {sample_name}")
        else:
            sample_name = df[0][ind]
            sample_unique = df[1][ind]
            if os.path.isfile(sample_unique):
                samples[sample_name] = [sample_unique]
                print("File for sample", sample_name, "identified.")
            else:
                print("Problem in reads file for sample", sample_name)
                print("Please specify correct file paths.")
                raise Exception(f"Problem in reads file for sample {sample_name}")
    return samples

def generate_config_file(samples, args):
    data_type = args["data_type"]
    run_name = args["run_name"]
    with open(args["config_file"], "w") as f:

        samples_string = "samples:\n"
        f.write(samples_string)
        for sample in samples.keys():
            if data_type == "illumina":
                my_sample = f"    sample-{str(sample)}: {samples[sample][0]} {samples[sample][1]}\n"
            else:
                my_sample = f"    sample-{str(sample)}: {samples[sample][0]}\n"
            f.write(my_sample)

        data_line = f"data: {data_type}\n"
        f.write(data_line)

        reference_line = f"reference: {args['reference']}\n"
        f.write(reference_line)

        scheme_line = f"scheme: {args['primer_scheme']}\n"
        f.write(scheme_line)

        minimum_depth_line = f"minimum_depth: {str(args['minimum_coverage'])}\n"
        f.write(minimum_depth_line)

        threads_line = f"threads: {str(args['threads'])}\n"
        f.write(threads_line)

        workflow_path_line = f"workflow_path: {sys.path[0]}\n"
        f.write(workflow_path_line)

        output_line = f"output: {args['output']}/{run_name}/\n"
        f.write(output_line)

        if args["data_type"] == "illumina":
            adapters_line = f"adapters: {args['adapters']}\n"
            f.write(adapters_line)

            minimum_length_line = (
                f"minimum_length: {str(args['minimum_read_length'])}\n"
            )
            f.write(minimum_length_line)

            trim_line = f"trim: {str(args['trim'])}\n"
            f.write(trim_line)

    print("The config file was generated. -> ", args["config_file"])

    return


def main(args):
    samples = validate_args(args)
    generate_config_file(samples, args)
    
    if args["create_config_only"]:
        print("Finished.")
        return 0

    config = args["config_file"]
    cores = args["threads_total"]
    target_rule = "all"

    thisdir = os.path.abspath(os.path.dirname(__file__))
    workflow_path = os.path.join(thisdir, 'scripts',f"consensus_{args['data_type']}.smk")

    successful = snakemake(
        workflow_path, configfiles=[config], cores=cores, targets=[target_rule]
    )

    return 0 if successful else 1
