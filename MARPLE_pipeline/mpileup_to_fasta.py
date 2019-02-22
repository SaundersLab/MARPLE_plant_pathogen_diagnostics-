#!/usr/bin/env python

import pandas as pd
import heapq
import sys

input_file = sys.argv[1]

with open(input_file) as snp_ratios_raw, open(input_file+".output","w") as output_file:
	for each_line in snp_ratios_raw:
		output_list = list()
		each_line_list = each_line.strip().split('\t')
		contig_id = each_line_list[0]
		position = each_line_list[1]
		ref_base = each_line_list[2]
		coverage = each_line_list[3]
		if len(each_line_list) > 4:
			alt_base = each_line_list[4]
			alt_base_freqs = each_line_list[5]
			alt_base_list = alt_base.split(',')
			alt_base_freqs_list = alt_base_freqs.split(',')
			alt_base_freqs_list = [float(element) for element in alt_base_freqs_list]
			sum_of_alt_base_freqs = str(float(sum(alt_base_freqs_list)))
			alt_bases = list(zip(alt_base_list, alt_base_freqs_list))
		else:
			alt_bases = 'NA'
			sum_of_alt_base_freqs = 'NA'
		output_list.extend((contig_id, position, coverage, ref_base))
		if type(alt_bases) == list and len(alt_bases) >= 1:
			output_list.append(str(alt_bases[0]))
		else:
			output_list.append('NA')
		if type(alt_bases) == list and len(alt_bases) >= 2:
			output_list.append(str(alt_bases[1]))
		else:
			output_list.append('NA')
		if type(alt_bases) == list and len(alt_bases) >= 3:
			output_list.append(str(alt_bases[2]))
		else:
			output_list.append('NA')
		if type(alt_bases) == list and len(alt_bases) >= 4 :
			output_list.append(str(alt_bases[3]))
		else:
			output_list.append('NA')
		output_list.append(sum_of_alt_base_freqs)
		output_file.write('\t'.join(output_list) + '\n')


coverage_cutoff = float(sys.argv[2])
min_allele_freq = float(sys.argv[3])
with open(input_file+".output") as snp_ratios_processed_file, open(input_file+".output.processed", 'w') as output_file_final:
	for each_line in snp_ratios_processed_file:
		each_line_list = each_line.strip().split('\t')
		output_file_final.write('\t'.join(each_line_list[:4]))
		if float(each_line_list[2]) < coverage_cutoff:
			output_file_final.write('\t' + 'NA' + '\t'+'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')
		else:
			if each_line_list[8] == 'NA':
				sum_allele_freq = 1
			else:
				sum_allele_freq = float(each_line_list[8])
			if len(each_line_list[4].strip('()').split(',')) == 2:
				base1, allelefreq1 = each_line_list[4].strip('()').split(',')
			else:
				base1 = None
				allelefreq1 = None
			if len(each_line_list[5].strip('()').split(',')) == 2:
				base2, allelefreq2 = each_line_list[5].strip('()').split(',')
			else:
				base2 = None
				allelefreq2 = None
			if len(each_line_list[6].strip('()').split(',')) == 2:
				base3, allelefreq3 = each_line_list[6].strip('()').split(',')
			else:
				base3 = None
				allelefreq3 = None
			if len(each_line_list[7].strip('()').split(',')) == 2:
				base4, allelefreq4 = each_line_list[7].strip('()').split(',')
			else:
				base4 = None
				allelefreq4 = None
			allelefreqs = [allelefreq1, allelefreq2, allelefreq3, allelefreq4]
			allelefreqs = [float(element) for element in allelefreqs if element != None]
			bases = [base1, base2, base3, base4]
			bases = [element for element in bases if element != 'NA']
			base_allelefreq = list(zip(bases,allelefreqs))
			base_allelefreq = sorted(base_allelefreq, key=lambda x: x[1], reverse=True)
			if len(allelefreqs) >= 1 and base_allelefreq[0][1] > (min_allele_freq * sum_allele_freq):
				output_file_final.write('\t' + base_allelefreq[0][0].strip('\'') +'\t' + str(base_allelefreq[0][1]))
			else:
				output_file_final.write('\t' + 'NA' + '\t' + 'NA')
			try:
				# if base_allelefreq[0][1] == base_allelefreq[1][1]:
				# 	print(base_allelefreq)
				if base_allelefreq[1][1]  > (min_allele_freq * sum_allele_freq):
					#BUG IN CODE - FIX!
					output_file_final.write('\t' + base_allelefreq[1][0].strip('\'') +'\t' + str(base_allelefreq[1][1]) +'\n')
				else:
					output_file_final.write('\t' + 'NA' + '\t' + 'NA' + '\n')
			except:
				output_file_final.write('\t' + 'NA' + '\t' + 'NA' + '\n')

with open(input_file+".output.processed") as input_file_for_getting_sequence, open(input_file+".output.processed.fasta",'w') as fasta_file:
	previous_sequence = None
	current_sequence = None
	for each_line in input_file_for_getting_sequence:
		each_line_list = each_line.strip().split('\t')
		current_sequence = each_line_list[0]
		if previous_sequence is None:
			fasta_file.write('>' + current_sequence + '\n')
			previous_sequence = current_sequence
		if current_sequence != previous_sequence:
			fasta_file.write('\n>' + current_sequence + '\n')
		base_set = set([each_line_list[4], each_line_list[6]])
		if each_line_list[4] == 'NA':
			fasta_file.write('?')
		elif each_line_list[4] != 'NA' and each_line_list[6] == 'NA':
			fasta_file.write(each_line_list[4])
		elif each_line_list[4] == each_line_list[6]:
			fasta_file.write(each_line_list[4])
		elif base_set == {'A', 'T'}:
			fasta_file.write('W')
		elif base_set == {'A', 'C'}:
			fasta_file.write('M')
		elif base_set == {'A', 'G'}:
			fasta_file.write('R')
		elif base_set == {'C', 'G'}:
			fasta_file.write('S')
		elif base_set == {'C', 'T'}:
			fasta_file.write('Y')
		elif base_set == {'G', 'T'}:
			fasta_file.write('K')
		previous_sequence = current_sequence
	fasta_file.write('\n')