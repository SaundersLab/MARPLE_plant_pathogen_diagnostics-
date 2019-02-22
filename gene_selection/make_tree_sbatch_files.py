#!/usr/bin/env python
import glob
list_of_alignment_files=glob.glob("filtered*.fasta")
for element in list_of_alignment_files:
    with open('slurm_submit_'+element.split('3000_')[1]+'.sh','w') as output_file:
        output_file.write('''#!/bin/sh
raxmlHPC-PTHREADS-SSE3 -T 10 -s {alignment} -m GTRGAMMA -n {tree_base} -p 100'''.format(alignment=element, tree_base=element.split('3000_')[1]+'.tree'))
