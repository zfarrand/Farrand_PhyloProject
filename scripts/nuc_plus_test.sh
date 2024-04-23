#!/bin/bash

SAMPLES="KG103
KG116
KG145
KG169
KG196
KG210
KG234
KG253
KG264
KG274
KG278
KG296
KG298
KG316
KG331
KG338
KG343
KG352
KG366
KG376
KG382
KG385
KG387
KG391
KG404
KG428
KG433
KG454
KG458
KG482
KG488
KG503
KG516
KG524
KG533
KG545
KG546
MSB285669
MSB285671
MSB285675
MSB285684"

for i in $SAMPLES; do
awk -F'\t' '{print $1}' coordinates/plus/XM_004578932.1_coord.txt | samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna -r - | bcftools consensus -s "$i" 41_norm_dp4_renamed.vcf.gz -I -M "N" > nuc_fasta/plus/XM_004578932.1/"$i"_XM-004578932.1.fa
done
#Then concatenate the fasta records for each exon within each fasta file
for infile in nuc_fasta/plus/XM_004578932.1/*.fa; do j=$(basename ${infile} .fa); sed -s '1h;/>/d;H;$!d;x;s/\n//2g' nuc_fasta/plus/XM_004578932.1/"$j".fa > nuc_fasta/plus/XM_004578932.1/"$j"_cat.fa; done
#Then did this to get the sample names into the fasta header
for infile in nuc_fasta/plus/XM_004578932.1/*_cat.fa; do f=$(basename ${infile} _cat.fa); sed -i "s/^>/>${f%_*}_/" nuc_fasta/plus/XM_004578932.1/"$f"_cat.fa; done
#Then combined each sample's fasta into multisample gene alignments. 
cat nuc_fasta/plus/XM_004578932.1/*XM-004578932.1_cat.fa > nucfasta/cat/XM_004578932.1.fa
#Then rename the sample headers in the fasta file so the genes all match and with more descriptive names
for infile in nuc_fasta/cat/*.fa; do g=$(basename ${infile} .fa); awk 'NR%2==0' nuc_fasta/cat/"$g".fa | paste -d'\n' nuc_fasta_rename.txt - > nuc_fasta/och_final/"$g"_renamed.fa; done