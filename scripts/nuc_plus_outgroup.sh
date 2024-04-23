#!/bin/bash

#Do a blast search for each pika gene to the rabbit genome, saving the rabbit sequence
for infile in coordinates/plus/*.txt; do base=$(basename ${infile} _coord.txt);
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna -r coordinates/plus/${base}_coord.txt | blastn -query - -db refgenome/OryCun2 -outfmt "6 sseq" -out OryCun2/${base}.fa 
done
#Then concatenate everything to the same line, remove spaces, and add a fasta header
for infile in OryCun2/*.fa; do base=$(basename ${infile} .fa);
echo $(cat OryCun2/${base}.fa) | sed 's/ //g' | sed -e '1i\>Oryctolagus_cuniculus' > OryCun2/final_nuc/${base}.fa
done
#Then add the rabbit sequence to the pikas fasta records
#Haven't checked this yet
for infile in nuc_fasta/final/*_renamed.fa > cat nuc_fasta/final/${base}_renamed*.fa OryCun2/final_nuc/${base}.fa > nuc_fasta/och_ory/${base}.fa; done

