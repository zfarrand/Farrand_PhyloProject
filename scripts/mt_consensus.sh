#!/bin/bash

FILES="KG103
KG116
KG145
KG169
KG196
KG210
KG228
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
MSB285674
MSB285680
MSB285682
MSB288917"

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
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:2745-3699 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/ND1/"$i"_ND1.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:3908-4949 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/ND2/"$i"_ND2.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:5333-6874 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/COX1/"$i"_COX1.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:7018-7701 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/COX2/"$i"_COX2.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:7773-7979 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/ATP8/"$i"_ATP8.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:7934-8613 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/ATP6/"$i"_ATP6.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:8614-9397 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/COX3/"$i"_COX3.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:9468-9813 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/ND3/"$i"_ND3.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:9883-10179 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/ND4L/"$i"_ND4L.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:10173-11550 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/ND4/"$i"_ND4.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:11749-13560 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/ND5/"$i"_ND5.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:13556-14080 -i | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/ND6/"$i"_ND6.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:14153-15290 | bcftools consensus -s "$i" 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz -M "N" > mt_fasta_dp3/CYTB/"$i"_CYTB.fa
done
#Then did this to get the sample names into the fasta header
for g in $GENES; do for infile in mt_fasta_dp3/"$g"/*.fa; do f=$(basename ${infile} .fa); sed -i "s/^>/>${f%_*}_/" mt_fasta_dp3/"$g"/"$f".fa; done; done
#Then combined each sample's fasta into multisample gene alignments. 
for g in $GENES; do cat mt_fasta_dp3/"$g"/*.fa > mt_fasta_dp3/"$g".fa; done
#Make fasta records into single line
for infile in mt_fasta_dp3/*.fa; do g=$(basename ${infile} .fa); awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' mt_fasta_dp3/"$g".fa > mt_fasta_dp3/"$g"_SL.fa; done
#Then rename the sample headers in the fasta file so the genes all match and with more descriptive names
for infile in mt_fasta_dp3/*_SL.fa; do g=$(basename ${infile} _SL.fa); awk 'NR%2==0' mt_fasta_dp3/"$g"_SL.fa | paste -d'\n' mt_fasta_rename.txt - > mt_fasta_dp3/final/"$g"_renamed.fa; done