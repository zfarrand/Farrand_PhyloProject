#!/bin/bash

TRANSCRIPTS="XM_004594296.2
XM_004599960.3
XM_004600005.1
XM_004600054.3
XM_012927943.1
XM_036492244.1
XM_036492781.1
XM_036497236.1"

for i in $TRANSCRIPTS
do
grep "$i" GCF_014633375.1_OchPri4.0_genomic.gff | awk -F'\t' -v OFS="\t" 'NR>1{ split($9,a,";");print $1,$4,$5,$7,a[6]}'| grep product | awk '{ print $1 ":" $2 "-" $3 "\t" $4 "\t" $5 }' > coordinates/minus/"$i"_coord.txt
done