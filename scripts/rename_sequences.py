#!/usr/bin/env python

import re

# Maybe include options to parse metadata
# and include in sequence names from here.

input = snakemake.input[0]
output = snakemake.output[0]
#print(input)

sample_name = re.sub(".+/","",input)
sample_name = re.sub(".consensus.fasta","",sample_name)
#print(sample_name)

with open(input) as f:
    lines = f.readlines()

sequence = ""

for line in lines:
    if line.startswith(">"):
        header = line
    else:
        sequence = sequence + line

with open(output, 'w') as f:
    renamed_sequence = ">" + sample_name + "\n" + sequence
    f.write(renamed_sequence)

exit()