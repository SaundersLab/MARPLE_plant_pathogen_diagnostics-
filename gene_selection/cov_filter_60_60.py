#!/usr/bin/env python

 
#import statements
from __future__ import division
from Bio import SeqIO
from collections import defaultdict, Counter

#function definitions
def coverage_filter(list_of_fasta_files, sequence_coverage_cutoff, sample_coverage_cutoff):
    #Get a way of opening all files together and reading one sequence at a time
    fasta_file_dictionary = { fasta_file:SeqIO.index(fasta_file,"fasta") for fasta_file in list_of_fasta_files }
    representative_fasta_file = list_of_fasta_files[0]
    list_of_sequences = sorted([ key for key in fasta_file_dictionary[representative_fasta_file].keys() ])
    coverage_dictionary = defaultdict(list)
    list_of_sequences_to_keep = []
    for each_sample in list_of_fasta_files:
        for each_sequence in list_of_sequences:
            if 1-(fasta_file_dictionary[each_sample][each_sequence].seq.count('N')/len(str(fasta_file_dictionary[each_sample][each_sequence].seq))) > sequence_coverage_cutoff:
                coverage_dictionary[each_sample].append(each_sequence)
    for each_sequence in list_of_sequences:
        number_of_samples_containing_this_sequence = len([ sample_name for sample_name in list_of_fasta_files if each_sequence in coverage_dictionary[sample_name]])
        if number_of_samples_containing_this_sequence/len(list_of_fasta_files) > sample_coverage_cutoff:
            list_of_sequences_to_keep.append(each_sequence)
    for each_sample in list_of_fasta_files:
        with open(each_sample+'.60_60_coverage_filtered.fa','w') as cov_filtered_file:
            for record in list_of_sequences_to_keep:
                cov_filtered_file.write(fasta_file_dictionary[each_sample][record].format("fasta"))



import glob

list1 = glob.glob("*cds.fa")

coverage_filter(list1,0.60,0.60)
