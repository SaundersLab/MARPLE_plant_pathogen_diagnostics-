#!/usr/bin/env python
import pandas as pd
input_table = pd.read_table('./gene_site_info.out.third_codon.list.sorted.tsv', header = 0)
lengths_table = pd.read_table('./gene_lengths.list',header = 0)
counts_table=input_table.groupby('Gene', as_index=False).count().sort_values("Position", ascending=False)
snps_per_length_table = counts_table.merge(lengths_table, on='Gene')
snps_per_length_table['snps_per_base']=snps_per_length_table['Position']/(snps_per_length_table['Length']/3.0)
sorted_snps_per_length_table = snps_per_length_table.sort_values('snps_per_base', ascending=False)
filtered_snps_table_1000_3000 = sorted_snps_per_length_table[(sorted_snps_per_length_table['Length']>1000) & (sorted_snps_per_length_table['Length']<3000)]
filtered_snps_table_1000_3000[filtered_snps_table_1000_3000['snps_per_base'] >= 0.005]
range_list = range(10, 1205, 5)
new_range = [ i / 10000.0 for i in range_list ]
#list_of_snp_values = [0.001]
list_of_snp_values = new_range
for each_cutoff in list_of_snp_values:
    filtered_snps_table_1000_3000[filtered_snps_table_1000_3000['snps_per_base']>=each_cutoff].to_csv('filtered_snps_length_1000to3000.snps_per_base_value_of_'+str(each_cutoff)+'.list', index=False, sep='\t')
    new_filtered_table = filtered_snps_table_1000_3000[filtered_snps_table_1000_3000['snps_per_base']>=each_cutoff]
    list_of_filtered_genes = list(new_filtered_table['Gene'])
    original_table_filtered = input_table[input_table['Gene'].isin(list_of_filtered_genes)]
    original_table_filtered.to_csv('filtered_table_length_1000_3000_snp_value_'+str(each_cutoff)+'.list.forfasta', index=False, sep='\t')
