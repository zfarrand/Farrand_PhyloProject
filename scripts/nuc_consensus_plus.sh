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
MSB285676
MSB285682
MSB285684
MSB293477
MSB293654"

for i in $SAMPLES; do for infile in coordinates/plus/*.txt; do base=$(basename ${infile} _coord.txt);
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna -r coordinates/plus/${base}_coord.txt | bcftools consensus -s "$i" 45_dp3.recode.vcf.gz -I -M "N" > nuc_fasta/plus/${base}/${base}_"$i".fa
done; done
#Then did this to get the sample names into the fasta header
for infile in coordinates/plus/*.txt; do base=$(basename ${infile} _coord.txt); for f in nuc_fasta/plus/${base}/*.fa; do sed -i "s/^>/>${f%_*}_/" nuc_fasta/plus/${base}/"$f"; done; done
#Then combined each sample's fasta into multisample gene alignments. 
for infile in coordinates/plus/*.txt; do base=$(basename ${infile} _coord.txt); cat ${base}*.fa > ${base}.fa; done