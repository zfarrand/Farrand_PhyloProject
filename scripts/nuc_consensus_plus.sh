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

for i in $SAMPLES; do for infile in coordinates/plus/*.txt; do f=$(basename ${infile} _coord.txt); awk -F'\t' '{print $1}' coordinates/plus/"$f"_coord.txt | samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna -r - | bcftools consensus -s "$i" 41_norm_dp4_renamed_ig5.vcf.gz -I -M "N" > nuc_fasta/plus/"$f"/"$i"_"$f".fa
done; done
#Then concatenate the fasta records for each exon within each fasta file
for infile in coordinates/plus/*.txt; do f=$(basename ${infile} _coord.txt); for infile in nuc_fasta/plus/"$f"/*.fa; do j=$(basename ${infile} .fa); sed -s '1h;/>/d;H;$!d;x;s/\n//2g' nuc_fasta/plus/"$f"/"$j".fa > nuc_fasta/plus/"$f"/"$j"_cat.fa
done; done
#Then did this to get the sample names into the fasta header
for infile in coordinates/plus/*.txt; do f=$(basename ${infile} _coord.txt); for infile in nuc_fasta/plus/"$f"/*_cat.fa; do j=$(basename ${infile} _cat.fa); sed -i "s/^>/>${f%_*}_/" nuc_fasta/plus/"$f"/"$j"_cat.fa
done; done
#Then combined each sample's fasta into multisample gene alignments. 
for infile in coordinates/plus/*.txt; do f=$(basename ${infile} _coord.txt); cat nuc_fasta/plus/"$f"/*_cat.fa > nuc_fasta/cat/"$f".fa
done
#Then rename the sample headers in the fasta file so the genes all match and with more descriptive names
for infile in coordinates/plus/*.txt; do f=$(basename ${infile} _coord.txt); awk 'NR%2==0' nuc_fasta/cat/"$f".fa | paste -d'\n' nuc_fasta_rename.txt - > nuc_fasta/och_final/"$f"_renamed.fa; done