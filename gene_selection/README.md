### Prerequisites
``` 
python-3.6.3
biopython-1.70
pandas-0.23.3
RAxML-8.0.12
```
### Protocol

Consensus CDS sequences for each RNAseq/genomic dataset were obtained as previously described (Bueno-Sancho et al. 2017) and placed in the folder where you will run these scripts

These consensus sequence files were then filtered for coverage of greater than 60% in each sample and in each gene
```
./cov_filter_60_60.py
```

The polymorphic sites were obtained using the following script
```
./gene_site_info_in_memory_filtered.py *.60_60_coverage_filtered.fa
```
Get third codon sites
```
./get_3rd_codon_sites.py gene_site_info.out
```
Sort snp list
```
./sort_snp_list.py gene_site_info.out.third_codon.list
```
Get lengths of cds sequences for your genome of interest
```
./write_gene_lengths.py PST130_cds.fa
```
Get filtered table
```
./from_sorted_snp_list_to_filtered_table.py
```
Get alignment from filtered table
```
for f in filtered_table_length_1000_3000_snp_value_0.*; do ./whole_genes_write_alignment_from_SNP_list.py $f *.60_60_coverage_filtered.fa; done
```
RAxML was then used to make trees using the alignments for the different SNP values
```
./make_tree_sbatch_files.py
```
Trees for the different SNP cutoffs were then annotated and visualized using dendroscope and MEGA
