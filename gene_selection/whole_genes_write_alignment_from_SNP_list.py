#!/usr/bin/env python

from Bio import SeqIO
import sys


with open(sys.argv[1]) as list_file:
	lines_in_list_file = list_file.readlines()
	set_of_lines_in_list_file = list(set([element.split('\t')[0] for element in lines_in_list_file]))
list_of_fasta_files = sys.argv[2:]
current_fasta_header = None
fasta_file_dictionary = { fasta_file:SeqIO.to_dict(SeqIO.parse(fasta_file,"fasta")) for fasta_file in list_of_fasta_files }
with open('alignment_from_'+sys.argv[1]+'.fasta','w') as output_fasta_file, open('third_codon_alignment_from_'+sys.argv[1]+'.fasta','w') as third_codon_output_fasta_file:
	for fasta_file in list_of_fasta_files:
		if fasta_file == list_of_fasta_files[0]:
			output_fasta_file.write('>'+fasta_file+'\n')
			third_codon_output_fasta_file.write('>'+fasta_file+'\n')
		else:
			output_fasta_file.write('\n>'+fasta_file+'\n')
			third_codon_output_fasta_file.write('\n>'+fasta_file+'\n')
		for each_line in set_of_lines_in_list_file:
			line_header = each_line
			if line_header.startswith('Gene'):
				continue
			else:
				output_fasta_file.write(str(fasta_file_dictionary[fasta_file][line_header].seq))
				third_codon_output_fasta_file.write(str(fasta_file_dictionary[fasta_file][line_header].seq[2::3]))
#				output_fasta_file.write(str(line_position)+str(fasta_file_dictionary[fasta_file][line_header].seq)[line_position-1])
