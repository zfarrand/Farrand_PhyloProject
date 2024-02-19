#!/bin/bash

SAMPLES="MSB285669
MSB293477
MSB293654
MSB285671
MSB285676
MSB285682
MSB285675
MSB285684
KG376
KG378
KG454
KG458
KG338
KG331
KG316
KG488
KG533
KG545
KG524
KG274
KG404
KG428
KG298
KG433
KG264
KG482
KG296
KG228
KG278
KG196
KG253
KG391
KG169
KG210
KG387
KG343
KG503
KG385
KG382
KG546
KG352
KG516
KG145
KG366
KG116
KG103"

for i in $SAMPLES; do; for infile in coordinates/plus/*.txt; do; base=$(basename ${infile} .txt);
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna -r ${base}.txt | bcftools consensus -s "$i" nuc_46.recode.vcf.gz -I -M "N" > nuc_fasta/${base}_"$i".fa
done
#Then did this to get the sample names into the fasta header
for f in nuc_fasta/*.fa; do; sed -i "s/^>/>${f%_*}_/" "$f"; done
#Then combined each sample's fasta into multisample gene alignments. 
for g in GENES; do; cat "$g"*.fa > "$g".fa; done
#Then align with MAFFT 
for g in GENES; do; mafft "$g".fa > "$g"/"$g".msa; done