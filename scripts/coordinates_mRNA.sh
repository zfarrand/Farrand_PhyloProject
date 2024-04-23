#!/bin/bash

TRANSCRIPTS="XM_004599994.1
XM_012928225.1
XM_012930006.1
XM_036492333.1
XM_036492437.1
XM_036493339.1
XM_036495691.1
XM_036495813.1
XM_004586853.1"

for i in $TRANSCRIPTS
do
grep "$i" GCF_014633375.1_OchPri4.0_genomic.gff | awk -F'\t' -v OFS="\t" 'NR>1{ split($9,a,";");print $1,$4,$5,$7,a[6]}'| grep mRNA | awk '{ print $1 ":" $2 "-" $3 "\t" $4 "\t" $5 }' > coordinates/minus/"$i"_coord.txt
done