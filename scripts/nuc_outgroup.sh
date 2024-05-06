#!/bin/bash

#Make fasta files for the pika genes
for infile in coordinates/all/*_coord.txt; do f=$(basename ${infile} _coord.txt); awk -F'\t' '{print $1}' coordinates/all/"$f"_coord.txt | samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna -r - > coordinates/gene_fasta/"$f"_ochpri4.fa; done

#Concatenate these exons into a single-line fasta
for infile in coordinates/gene_fasta/*_ochpri4.fa; do f=$(basename ${infile} _ochpri4.fa); sed -s '1h;/>/d;H;$!d;x;s/\n//2g' coordinates/gene_fasta/"$f"_ochpri4.fa > coordinates/gene_fasta/"$f"_ochpri4_cat.fa; done

#Run a blast search with these genes on the rabbit genome
for infile in coordinates/gene_fasta/*_ochpri4_cat.fa; do f=$(basename ${infile} _ochpri4_cat.fa); blastn -query coordinates/gene_fasta/"$f"_ochpri4_cat.fa -db refgenome/Orycun2 -outfmt "6 sseq" -max_target_seqs 1 -task blastn -out OryCun2/"$f".fa 
done

#Then concatenate everything to the same line, remove spaces, and add a fasta header
for infile in OryCun2/*.fa; do j=$(basename ${infile} .fa);
echo $(cat OryCun2/"$j".fa) | sed 's/ //g' | sed -e '1i\>Oryctolagus_cuniculus' > OryCun2/final_nuc/"$j".fa
done

#Then add the rabbit sequence to the pikas fasta records
for infile in nuc_fasta/och_final/*_renamed.fa; do f=$(basename ${infile} _renamed.fa); cat nuc_fasta/och_final/"$f"_renamed.fa OryCun2/final_nuc/"$f".fa > nuc_fasta/och_ory/"$f".fa; done

############################

awk -F'\t' '{print $1}' coordinates/all/NM_001310092.1_coord.txt | samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna -r - | blastn -query - -db refgenome/Orycun2 -task blastn -outfmt "6 sseq" -out outgroup_troubleshooting/NM_001310092.1.fa 

#Maybe I need to put the samtools faidx region into a single line file before blasting? 

awk -F'\t' '{print $1}' coordinates/all/NM_001310092.1_coord.txt | samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna -r - > outgroup_troubleshooting/NM_001310092.1_ochpri4.fa 

sed -s '1h;/>/d;H;$!d;x;s/\n//2g' outgroup_troubleshooting/NM_001310092.1_ochpri4.fa > outgroup_troubleshooting/NM_001310092.1_ochpri4_cat.fa 

blastn -query outgroup_troubleshooting/NM_001310092.1_ochpri4_cat.fa -db refgenome/Orycun2 -outfmt "6 sseq" -out outgroup_troubleshooting/NM_001310092.1_ory.fa 
