#!/bin/bash

#arguments
fastafile=$2
fastqfile=$1
gff_file=final_poolD_genes.gff3
samfile=$fastqfile.$fastafile.sam
bamfile=$samfile.bam
SET='\033[0m'
echo -e "[1m[4mPROCESSING SAMPLE - $1 $SET"


#Porechop for filtering and trimming the reads
porechop -i $fastqfile -o $fastqfile.porechopped -t 8 > $fastqfile.porechop.log

head -5 $fastqfile.porechop.log

mv $fastqfile $fastqfile.initial
mv $fastqfile.porechopped $fastqfile
echo -e "[1m[4mFILTERING AND TRIMMING IS DONE$SET"

#Mapping
bwa index $fastafile
bwa mem $fastafile $fastqfile > $samfile
#minimap2 -x map-ont -a -L $fastafile $fastqfile > $samfile
samtools sort -O BAM -o $bamfile $samfile
samtools index $bamfile

echo -e "[1m[4mMAPPING IS DONE$SET"

#SNP calling
samtools mpileup -BQ 0 -x -d 200000 -aa -f $fastafile $bamfile > $bamfile.mpileup
cat $bamfile.mpileup | python compsnps_pipe1_sampileup.py > $bamfile.mpileup.snp_ratios.txt

echo -e "[1m[4mSNP CALLING IS DONE$SET"

#Get fasta file from SNP called mpileup files
#Usage: ./mpileup_to_fasta.py SNP_ratios_file Coverage_cutoff Allele_freq_cutoff
python mpileup_to_fasta.py $bamfile.mpileup.snp_ratios.txt 20 0.25

echo -e "[1m[4mPRELIMINARY FASTA IS READY$SET"

mv $bamfile.mpileup.snp_ratios.txt.output.processed.fasta $bamfile.mpileup.snp_ratios.txt.output.processed.fasta.genes

#GFF filter step to get cds from genes
gffread $gff_file -g $bamfile.mpileup.snp_ratios.txt.output.processed.fasta.genes -w $bamfile.mpileup.snp_ratios.txt.output.processed.fasta

#Sort fasta file
seqkit sort -n $bamfile.mpileup.snp_ratios.txt.output.processed.fasta > $bamfile.mpileup.snp_ratios.txt.output.processed.fasta.sorted

#Concatenate fasta file
union -sequence $bamfile.mpileup.snp_ratios.txt.output.processed.fasta.sorted -outseq $bamfile.mpileup.snp_ratios.txt.output.processed.fasta.sorted.concatenated

#Rename sequence header with filename
awk '/^>/ {gsub(/.fa(sta)?$/,"",FILENAME);printf(">%s\n",FILENAME);next;} {print}' $bamfile.mpileup.snp_ratios.txt.output.processed.fasta.sorted.concatenated > $bamfile.mpileup.snp_ratios.txt.output.processed.fasta.sorted.concatenated.renamed.fasta

echo -e "[1m[4mFASTA CONCATENATION IS DONE$SET"


#move files to the appropriate folders
mkdir -p intermediate
mkdir -p final_output
mv $bamfile.mpileup.snp_ratios.txt.output.processed.fasta.sorted.concatenated.renamed.fasta final_output/
mv $bamfile $bamfile.* intermediate

echo -e "[1m[4mFINISHED PROCESSING $2
$SET"

echo -e "[1m[4mTHE FINAL OUTPUT IS IN THE "final_output" folder for this particular sample - $2
$SET"

echo -e "[1m[4mNow you can move to the "final_output" folder and copy the file containing sequences from reference isolates that you want in your tree (name this file as 'reference_for_tree.fasta') and then run tree.sh to start constructing the tree
$SET"

