#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
from collections import defaultdict
import sys

def write_sequence_information(list_of_fasta_files, output_file_prefix):
	fasta_file_dictionary = { fasta_file:SeqIO.to_dict(SeqIO.parse(fasta_file,"fasta")) for fasta_file in list_of_fasta_files }
	representative_fasta_file = list_of_fasta_files[0]
	list_of_sequences = sorted([ key for key in fasta_file_dictionary[representative_fasta_file].keys() ])
	total_seq_dict = defaultdict()
	for each_sequence in list_of_sequences:
		each_seq_dict = defaultdict(list)
		final_seq_dict = defaultdict(list)
		unique_bases_in_position = []
		for position, dead_base in enumerate(str(fasta_file_dictionary[list_of_fasta_files[0]][each_sequence].seq)):
			for each_sample in list_of_fasta_files:
				each_seq_dict[str(position+1)].append(str(fasta_file_dictionary[each_sample][each_sequence].seq)[position])
			if each_seq_dict[str(position+1)].count('N') < (0.2*len(each_seq_dict[str(position+1)])):
				unique_bases_in_position = set(each_seq_dict[str(position+1)])
			if unique_bases_in_position != []:
				if 'N' in unique_bases_in_position: unique_bases_in_position.remove('N')
			if not(unique_bases_in_position is None):
				if len(unique_bases_in_position) > 1:
					final_seq_dict[str(position+1)] = list(unique_bases_in_position)
		each_seq_dict = None
		total_seq_dict[each_sequence]=final_seq_dict
		final_seq_dict = None
	with open(output_file_prefix, 'w') as output_table:
		output_table.write('Gene\tPosition\tPossible_bases\n')
		for sequence in total_seq_dict:
			for base_position in total_seq_dict[sequence]:
				output_table.write(sequence+'\t'+base_position+'\t'+','.join(total_seq_dict[sequence][base_position])+'\n')

if __name__ == '__main__':
	write_sequence_information(sys.argv[1:],'gene_site_info.out')
