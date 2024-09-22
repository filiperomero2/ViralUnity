#!/usr/bin/env python

"""
This scripts dynamically creates config files and run 
the viralunity consensus snakemake pipeline.
Filipe Moreira - 2024/09/21
"""

import os
import sys
import argparse
import pandas as pd
import datetime
from snakemake import snakemake

def get_args():
    
    parser = argparse.ArgumentParser(
    description='A script to generate config files and run the viralunity consensus pipeline',
    usage='''viralunity_consensus.py [args]''')

    parser.add_argument('--data-type', type = str,
    help='Sequencing data type (illumina or nanopore).',
    choices = ['illumina','nanopore'],
    nargs='?',required = True)

    parser.add_argument('--sample-sheet',
    help='Complete path for a csv file with samples data paths and metadata',
    required = True)
 
    parser.add_argument('--config-file', 
    help='Complete path for input (viralunity config) file to be created.',
    required = True)
    
    parser.add_argument('--output', 
    help='Complete path for output directory to be created by viralunity.',
    required = True)

    parser.add_argument('--run-name', type = str, 
    help='Name for the sequencing run (optional).',
    nargs='?',const=1, default='undefined')
    
    parser.add_argument('--reference', 
    help='Complete path for the reference genome in fasta format',
    required = True)

    parser.add_argument('--minimum-coverage',type = int,
    help='Minimum sequencing coverage for including base in consensus sequence (Default = 20)',
    nargs='?',const=1, default=20)

    parser.add_argument('--adapters', 
    help='Complete path for adapters sequences in fasta format [Illumina QC]')

    parser.add_argument('--minimum-read-length',type = int,
    help='Minimum read length threshold (Default = 50) [Illumina QC]',
    nargs='?',const=1, default=50)

    parser.add_argument('--trim',type = int,
    help="Number of bases to trim from the 5' end of reads (Default = 0) [Illumina QC]",
    nargs='?',const=1, default=0)

    parser.add_argument("--create-config-only",
    help="Only create config file, not running the workflow (boolean)",
    action='store_true')

    parser.add_argument('--threads',type = int,
    help='Number of available threads for individual tasks (Default = 1)',
    nargs='?',const=1, default=1)

    parser.add_argument('--threads-total',type = int,
    help='Number of available threads for the entire workflow (Default = 1)',
    nargs='?',const=1, default=1)

    args = vars(parser.parse_args())
    return args

def validate_args(args):
    
    if os.path.isfile(args['sample_sheet']):
        sample_sheet = args['sample_sheet']
        samples = validate_sample_sheet(sample_sheet,args)
    else:
        print("Sample sheet file does not exist, please verify.")
        exit()
    
    if(os.path.isfile(args['config_file'])):
        print("Config file already exists. Please specify a new one.")
        exit()

    if(os.path.isdir(args['output'])):
        print("Output directory already exists. Please specify a new one.")
        exit()

    if(os.path.isfile(args['reference'])):
        print("Reference sequence file exists.")
    else:
        print("Reference sequence file does not exist, please verify.")
        exit()

    if args['data_type'] == 'illumina':
        if(os.path.isfile(args['adapters'])):
            print("Illumina adapter sequences file exists.")
        else:
            print("Illumina adapter sequences file not found, please verify.")
            exit()

    print("All arguments were succesfully verified.")

    return(samples)

def validate_sample_sheet(sample_sheet,args): 
    df =  pd.read_csv(sample_sheet,header=None)
    samples = {}
    for ind in df.index:
        if args['data_type'] == 'illumina':
            # This assumes all illumina data are paired-end reads
            sample_name = df[0][ind]
            sample_R1 = df[1][ind]
            sample_R2 = df[2][ind]
            if os.path.isfile(sample_R1) and os.path.isfile(sample_R2):
                samples[sample_name] = [sample_R1,sample_R2]
                print("Files for sample",sample_name,"identified.")
            else:
                print("Problem in reads files for sample", sample_name)
                print("Please specify correct file paths.")
                exit()
        else:
            sample_name = df[0][ind]
            sample_unique = df[1][ind]
            if os.path.isfile(sample_unique):
                samples[sample_name] = [sample_unique]
                print("File for sample",sample_name,"identified.")
            else:
                print("Problem in reads file for sample", sample_name)
                print("Please specify correct file paths.")
                exit()
    return samples

def define_job_id(args):
    current_time = str(datetime.datetime.now()).split(' ')[0]
    id = args['run_name']
    job_id = 'job_' + str(current_time) + "_" + str(id)
    return job_id

def generate_config_file(samples,args):
    with open(args['config_file'], 'w') as f:

        samples_string = 'samples:' + "\n"
        f.write(samples_string)
        for sample in samples.keys():
            if args['data_type'] == 'illumina':
                my_sample = "    sample-" + str(sample) + ": " + samples[sample][0] + " " + samples[sample][1] + "\n"
            else:
                my_sample = "    sample-" + str(sample) + ": " + samples[sample][0] + "\n"
            f.write(my_sample)

        data = "data: " + args['data_type'] + "\n"
        f.write(data)

        reference = "reference: " + args['reference'] + "\n"
        f.write(reference)

        minimum_depth = "minimum_depth: " + str(args['minimum_coverage']) + "\n"
        f.write(minimum_depth)

        threads = "threads: " + str(args['threads']) + "\n"
        f.write(threads)

        workflow_path = "workflow_path: " + sys.path[0] + "\n"
        f.write(workflow_path)

        job_id = define_job_id(args)
        output = "output: " + args['output'] + "/" + job_id + "/" + "\n"
        f.write(output)

        if args['data_type'] == 'illumina':
            adapters = "adapters: " + args['adapters'] + "\n"
            f.write(adapters)

            minimum_length = "minimum_length: " + str(args['minimum_read_length']) + "\n"
            f.write(minimum_length)

            trim = "trim: " + str(args['trim']) + "\n"
            f.write(trim)
        
    print("The config file was generated. -> ", args['config_file'])

    return

def main():
    args = get_args()
    samples = validate_args(args)
    generate_config_file(samples,args)

    config = args['config_file']
    cores = args['threads_total']
    target_rule = "all"

    if args['data_type'] == 'illumina':
        workflow_path = sys.path[0] + '/consensus_illumina.smk'
    else:
        workflow_path = sys.path[0] + '/consensus_nanopore.smk'
    
    if args['create_config_only']:
        print('Finished.')
        status = True
    else:
        status = snakemake(workflow_path, configfiles=[config], cores=cores,targets=[target_rule])

    if status:
        return 0
    else:
        return 1

if __name__  == '__main__':
    sys.exit(main())

