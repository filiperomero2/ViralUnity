#!/usr/bin/env python

import os
import argparse
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def get_args():
    
    parser = argparse.ArgumentParser(
    description='A script to generate report with all mutations found in genome sequences.',
    usage='''mutation_mapper.py [args]''')
 
    parser.add_argument('--input', 
    help='Complete path for the input alignment in fasta file. It assumes the first sequence is the reference genome.',
    required = True)

    parser.add_argument('--output', 
    help='Output annotation file.',
    required = True)
    
    parser.add_argument('--annotation-file',
    help='Complete path for annotation (gff) file.',
    required = True)

    parser.add_argument('--flag',
    help='Flag used to retrieve geve positions [CDS,gene]',
    choices = ['CDS','gene'],
    nargs='?',const=1, default='gene')

    args = vars(parser.parse_args())
    return args

def validate_args(args):

    if os.path.isfile(args['input']):
        print(f"Input alignment file identified -> {args['input']}")
    else:
        print(f"Input alignment file not identified, please verify.")
        exit()

    if os.path.isfile(args['output']):
        print(f"Output annotation report file already exists, please verify.")
        exit()
    
    if os.path.isfile(args['annotation_file']):
        print(f"Annotation input file identified -> {args['annotation_file']}")
    else:
        print(f"Annotation input file not identified, please verify.")
        exit()

    return

def get_annotation_info(args):
    df = pd.read_csv(args['annotation_file'],header=None,comment="#",sep="\t")
    
    genes = []
    metadata = tuple(df[df[2] == 'gene'][8])
    for info in metadata:
        my_match = re.search(r'Name=.*?;',info).group().split('=')[1].rstrip(';')
        genes.append(my_match)
    
    genes = tuple(genes)
    start = tuple(df[df[2] == 'gene'][3])
    end = tuple(df[df[2] == 'gene'][4])
    
    annotation_info = {}
    for i in range(0,len(genes)):
        annotation_info[genes[i]] = (start[i],end[i])

    return annotation_info

def read_and_annotate_sequences(args,annotation_info):
    
    data = {}
    genomic_mismatches = {}
    protein_mismatches = {}
    my_gene_data = {}

    with open(args['input']) as aln:
        for record in SeqIO.parse(aln, "fasta"):
            data[record.id] = record.seq
    
    ref_seq_name = list(data.keys())[0]
    ref_genome_seq = data[ref_seq_name]
    
    for seq_name in data.keys():
        if seq_name != ref_seq_name:
            
            my_genome_seq = data[seq_name]
            genomic_mismatches[seq_name] = compare_sequences(ref_genome_seq,my_genome_seq)
            
            my_gene_data[seq_name] = {}

            for gene in annotation_info.keys():

                # Eventually include an option to include a flag marking early stop codons
                my_gene_ref_seq = ref_genome_seq[annotation_info[gene][0]-1:annotation_info[gene][1]]
                my_gene_ref_seq = Seq(str(my_gene_ref_seq).replace('-','N')).translate()
                
                my_gene_seq = my_genome_seq[annotation_info[gene][0]-1:annotation_info[gene][1]]
                my_gene_data[seq_name][gene] = my_gene_seq
                my_gene_seq = Seq(str(my_gene_seq).replace('-','N')).translate()

                my_aa_mismatches = compare_sequences(my_gene_ref_seq,my_gene_seq)
                
                if len(my_aa_mismatches) == 0:
                    continue
                
                my_aa_mismatches = [gene + ':' + m  for m in my_aa_mismatches]
                
                if seq_name in protein_mismatches:
                    protein_mismatches[seq_name].append(my_aa_mismatches)
                else:
                    protein_mismatches[seq_name] = [my_aa_mismatches]
            
            protein_mismatches[seq_name] = sum(protein_mismatches[seq_name],[])
    
    path = args['output']
    path = "/".join(path.split('/')[:-1])

    if os.path.isdir(path) is False:
        os.mkdir(path)

    with open(args['output'],'w') as out:
        for seq_name in genomic_mismatches.keys():
            my_genomic_mismatches = str(genomic_mismatches[seq_name]).lstrip('[').rstrip(']').replace(' ','')
            my_protein_mismatches = str(protein_mismatches[seq_name]).lstrip('[').rstrip(']').replace(' ','')
            my_line_output_report = seq_name + "; " + my_genomic_mismatches + "; " + my_protein_mismatches + "\n"
            out.write(my_line_output_report)
    
    for gene in my_gene_data[seq_name].keys():
        out_gene = f"{path}/{gene}.fasta"
        with open(out_gene,'w') as out:
            for seq_name in my_gene_data.keys():  
                my_out = ">" + seq_name + "_" + gene + "\n" + str(my_gene_data[seq_name][gene]) + "\n"
                out.write(my_out)
    
def compare_sequences(ref_seq,my_seq):
    count = 0
    mutations = []
    for i in range(0,len(ref_seq)):
        ref_site = ref_seq[i]
        my_site = my_seq[i]
        # Indels not included in the report file
        if ref_site != my_site and my_site != 'N' and my_site != '-' and my_site != 'X' :
            count += 1
            mutation = ref_site + str(i) + my_site
            mutations.append(mutation)
    return(mutations)

def main():
    args = get_args()
    validate_args(args)
    annotation_info = get_annotation_info(args)
    read_and_annotate_sequences(args,annotation_info)
    print("The annotation output file was created. -> ", args['output'])

if __name__  == '__main__':
    main()
