#!/usr/bin/env python

"""
This scripts dynamically creates config files and run 
the viralunity metagenomics snakemake pipeline.
Filipe Moreira - 2024/02/14
"""

import os
import sys
import argparse
import pandas as pd
import datetime
from snakemake import snakemake

def get_args():
    
    parser = argparse.ArgumentParser(
    description='A script to generate config files and run the viralunity metagenomics pipeline',
    usage='''viralunity_meta.py [args]''')

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
    
    parser.add_argument('--kraken2-database', 
    help='Complete path for the kraken2 database directory',
    required = True)

    parser.add_argument('--krona-database', 
    help='Complete path for the krona taxonomic database',
    required = True)

    parser.add_argument('--adapters', 
    help='Complete path for Illumina adapters sequences in fasta format',
    required = True)

    parser.add_argument('--minimum-read-length',type = int,
    help='Minimum read length threshold (Default = 50)',
    nargs='?',const=1, default=50)

    parser.add_argument('--trim',type = int,
    help="Number of bases to trim from the 5' end of reads (Default = 0)",
    nargs='?',const=1, default=0)

    parser.add_argument("--remove-human-reads",
    help="Remove human reads from krona plot (boolean)",
    action='store_true')

    parser.add_argument("--remove-unclassified-reads",
    help="Remove unclassified reads from krona plot (boolean)",
    action='store_true')

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
        samples = validate_sample_sheet(sample_sheet)
    else:
        print("Sample sheet file does not exist, please verify.")
        exit()
    
    if(os.path.isfile(args['config_file'])):
        print("Config file already exists. Please specify a new one.")
        exit()

    if(os.path.isdir(args['output'])):
        print("Output directory already exists. Please specify a new one.")
        exit()

    if(os.path.isdir(args['kraken2_database'])):
        print("kraken2 database directory exists.")
    else:
        print("kraken2 database directory does not exist, please verify.")
        exit()

    if(os.path.isdir(args['krona_database'])):
        print("krona taxonomic database directory exists.")
    else:
        print("krona taxonomic database does not exist, please verify.")
        exit()

    if(os.path.isfile(args['adapters'])):
        print("Adapter sequences file exists.")
    else:
        print("Adapter sequence file does not exist, please verify.")
        exit()

    print("All arguments were succesfully verified.")

    return(samples)

def validate_sample_sheet(sample_sheet):
    df =  pd.read_csv(sample_sheet,header=None)
    samples = {}
    for ind in df.index:
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
            my_sample = "    sample-" + str(sample) + ": " + samples[sample][0] + " " + samples[sample][1] + "\n"
            f.write(my_sample)

        kraken2_database = "kraken2_database: " + args['kraken2_database'] + "\n"
        f.write(kraken2_database)

        krona_database = "krona_database: " + args['krona_database'] + "\n"
        f.write(krona_database)

        adapters = "adapters: " + args['adapters'] + "\n"
        f.write(adapters)

        minimum_length = "minimum_length: " + str(args['minimum_read_length']) + "\n"
        f.write(minimum_length)

        trim = "trim: " + str(args['trim']) + "\n"
        f.write(trim)

        if args['remove_human_reads']:
            remove_human_reads = "remove_human_reads: True" + "\n"
        else: 
            remove_human_reads = "remove_human_reads: False" + "\n"
        f.write(remove_human_reads)

        if args['remove_unclassified_reads']:
            remove_unclassified_reads = "remove_unclassified_reads: True" + "\n"
        else: 
            remove_unclassified_reads = "remove_unclassified_reads: False" + "\n"
        f.write(remove_unclassified_reads)

        threads = "threads: " + str(args['threads']) + "\n"
        f.write(threads)

        job_id = define_job_id(args)
        output = "output: " + args['output'] + "/" + job_id + "/" + "\n"
        f.write(output)

        print("The config file was generated. -> ", args['config_file'])

    return

def main():
    args = get_args()
    samples = validate_args(args)
    generate_config_file(samples,args)

    config = args['config_file']
    cores = args['threads_total']
    target_rule = "all"
    workflow_path = sys.path[0] + '/metagenomics.smk'

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
