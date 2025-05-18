#!/usr/bin/env python

# Create sample sheet from directory structure 
# Filipe Moreira - 2023/09/16

import os
import argparse
import glob
import sys

def get_args(args):
    parser = argparse.ArgumentParser(
    description='A script to generate a sample-sheet csv file from sequencing run directory structure',
    usage='''viralunity_create_samplesheet.py [args]''')
 
    parser.add_argument('--input', 
    help='Complete path for the directory containing sequencing data (0) or subdirectories with sequencing data (1)',
    required = True)
    
    parser.add_argument('--separator',
    help='Separator character to determine sample name.',
    choices = ['_','-','.'],
    nargs='?',const=1, default='-')

    parser.add_argument('--pattern',
    help='Pattern string unpaired reads file names.',
    choices = ['R1','barcode'],
    nargs='?',const=1, default='R1')

    parser.add_argument('--level', type = int,
    help='Directory level to perform sequencing data files search (0 or 1).',
    choices = [0,1],
    nargs='?',const=1, default=1)

    parser.add_argument('--output', 
    help='Output file name.',
    required = True)

    # Add separator and data type options
    # outdated script
    
    args = vars(parser.parse_args(args))
    
    return args

def validate_args(args):

    if os.path.isdir(args['input']) is False:
        print("Input directory does not exist, please verify.")
        raise Exception("Input directory does not exist.")
    else:
        print(f"Input data directory -> {args['input']}")

    if args['level'] == 1:
        print("Data will be searched for in its subdirectories.")
    else:
        print("Data will be searched for only in this directory.")
    
    if os.path.isfile(args['output']):
        print("Output file already exists, please verify.")
        raise Exception("Output file already exists.")
    
    return

def generate_sample_sheet(args):
    samples = {}
    if args['level'] == 1:
        directories = sorted(glob.glob(os.path.join(args['input'], '*')))
        for directory in directories:
            if os.path.isdir(directory):
                sample_name = directory.split('/')[-1].split(args['separator'])[0]
                sequencing_files_path = sorted(glob.glob(os.path.join(directory,'*')))
                if len(sequencing_files_path) == 2 and os.path.isfile(sequencing_files_path[0]) and os.path.isfile(sequencing_files_path[1]):
                    samples[sample_name] = sequencing_files_path
                    print(f"Two files for sample {sample_name} were identified. Proceeding...")
                elif len(sequencing_files_path) == 1 and os.path.isfile(sequencing_files_path[0]):
                    samples[sample_name] = sequencing_files_path
                    print(f"One file for sample {sample_name} was identified. Proceeding...")
                else:
                    print(f"Number of files found for sample {sample_name} is not one or two. Please verify.")
                    raise Exception(f"Number of files found for sample {sample_name} is not one or two. Please verify.")
    else:
        files = sorted(glob.glob(os.path.join(args['input'], f"*{args['pattern']}*")))
        for file in files:
            sample_name = file.split('/')[-1].split(args['separator'])[0]
            sequencing_files_path = sorted(glob.glob(os.path.join(args['input'], f"{sample_name}*")))
            if len(sequencing_files_path) == 2 and os.path.isfile(sequencing_files_path[0]) and os.path.isfile(sequencing_files_path[1]):
                samples[sample_name] = sequencing_files_path
                print(f"Two files for sample {sample_name} were identified. Proceeding...")
            elif len(sequencing_files_path) == 1 and os.path.isfile(sequencing_files_path[0]):
                    samples[sample_name] = sequencing_files_path
                    print(f"One file for sample {sample_name} was identified. Proceeding...")
            else:
                print(f"Number of files found for sample {sample_name} is not one or two. Please verify.")
                raise Exception(f"Number of files found for sample {sample_name} is not one or two. Please verify.")

    with open(args['output'], 'w') as f:
        for key in samples.keys():
            if len(samples[key]) == 2:
                sample_line = f"{key},{samples[key][0]},{samples[key][1]}\n"
                f.write(sample_line)
            else:
                sample_line = f"{key},{samples[key][0]}\n"
                f.write(sample_line)

def main():
    args = get_args(sys.argv[1:])
    validate_args(args)
    generate_sample_sheet(args)
    print("The sample sheet file was generated. -> ", args['output'])

if __name__  == '__main__':
    main()
