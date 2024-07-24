#!/usr/bin/env python

import os

raw_consensus = snakemake.input[0]
minimum_depth_table = snakemake.input[1]
output = snakemake.output[0]

if os.path.getsize(minimum_depth_table) == 0:
    print("Masking file is empty, copying raw consensus as final consensus.")
    copying_command = "cp " + raw_consensus + " " + output
    os.system(copying_command)
else:
    print("Low coverage sites identified, masking...")
    masking_command = "bedtools maskfasta -fi " + raw_consensus + " -fo " + output + " -bed " + minimum_depth_table
    os.system(masking_command)

exit()

