#!/bin/bash

FILES="MSB285669
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

GENES="ND1
ND2
COX1
COX2
ATP8
ATP6
COX3
ND3
ND4L
ND4
ND5
ND6
CYTB"

for i in $FILES
do
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:2745-3699 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/ND1/"$i"_ND1.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:3908-4949 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/ND2/"$i"_ND2.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:5333-6874 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/COX1/"$i"_COX1.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:7018-7701 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/COX2/"$i"_COX2.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:7773-7979 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/ATP8/"$i"_ATP8.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:7934-8613 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/ATP6/"$i"_ATP6.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:8614-9397 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/COX3/"$i"_COX3.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:9468-9813 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/ND3/"$i"_ND3.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:9883-10179 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/ND4L/"$i"_ND4L.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:10173-11550 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/ND4/"$i"_ND4.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:11749-13560 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/ND5/"$i"_ND5.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:13556-14080 -i | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/ND6/"$i"_ND6.fa
samtools faidx GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:14153-15290 | bcftools consensus -s "$i" mt_46_dp3.recode.vcf.gz -I -M "N" > mt_fasta_dp3/CYTB/"$i"_CYTB.fa
done
#Then did this to get the sample names into the fasta header
for g in GENES; do; for f in "$g"/*.fa; do; sed -i "s/^>/>${f%_*}_/" "$f"; done
#Then combined each sample's fasta into multisample gene alignments. 
for g in GENES; do; cat "$g"/*.fa > "$g".fa; done
#Then align with MAFFT 
for g in GENES; do; mafft "$g"/"$g".fa > "$g"/"$g".msa; done