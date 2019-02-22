#!/usr/bin/env python

import sys

with open(sys.argv[1]) as input_file, open(sys.argv[1]+'.third_codon.list','w') as output_file:
	lines = input_file.readlines()
	for eachline in lines:
		if eachline.startswith("Gene\t"):
			output_file.write(eachline)
		else:
			if int(eachline.split('\t')[1])%3==0:
				output_file.write(eachline)
			else:
				continue

