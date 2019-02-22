#!/usr/bin/env python

import pandas as pd
import sys

input_table = pd.read_table(sys.argv[1],header=0)
sorted_table=input_table.sort_values(['Gene', 'Position'])
sorted_table.to_csv(sys.argv[1]+'.sorted.tsv',index=False, sep='\t')

