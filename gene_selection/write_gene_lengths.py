#!/usr/bin/env python
from Bio import SeqIO
import sys

with open(sys.argv[1]) as input_file, open('gene_lengths.list','w') as gene_lengths_file:
    gene_lengths_file.write('Gene\tLength\n')
    for record in SeqIO.parse(input_file,'fasta'):
        gene_lengths_file.write(str(record.id)+'\t'+str(len(record.seq))+'\n')