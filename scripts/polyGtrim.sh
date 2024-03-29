#!/bin/bash

FILES="KG102_GGTGTACA-TAGCGTCT_L001
KG103_CCATCCGC-TACGCCTT_L001
KG108_CTTGCTAG-GGAGATGA_L001
KG111_TCGGATTC-TAGTTGCG_L001
KG114_TAGGAGCT-GTTAAGGC_L001
KG116_GTTGGCGT-GACTGCTG_L001
KG117_AACACCAC-CTCACACC_L001
KG128_CTGTACCA-AGTCTGTG_L001
KG136_GTCAGTCA-TAATGCCG_L001
KG139_GGTGTGAG-GTGGTGTT_L001
KG140_ATCTGACC-AAGTCGAG_L001
KG145_TCTACGCA-GGCTATTG_L001
KG148_CTTTCCCT-AGTAAGCC_L001
KG152_CCTATGCA-ATCGACTC_L001
KG153_TAATCTCG-GTGTGACA_L001
KG157_TACGGCAG-TTCATATC_L001
KG158_AGCGAGAT-GATACTGG_L001
KG169_TAGCTGGC-AGCCGGTA_L001
KG171_TGCCCATC-CCACACTT_L001
KG172_TTACCGAC-CGTATTCG_L001
KG191_CAGCAGTC-ACGCCAAC_L001
KG195_TCGTCTGA-TCAAGGAC_L001
KG196_TTCCAGGT-AAGCACTG_L001
KG201_GTCCGATC-GGATTGGC_L001
KG206_GTAAGGTG-GCAATTCG_L001
KG210_TATGACGT-GATTAGCG_L001
KG214_GAACGAAG-GGTGATTC_L001
KG215_ATTACCCA-CTACAGTG_L001
KG217_TCTAGGAG-AGGTCACT_L001
KG224_TGGTGAAG-GCTGTAAG_L001
KG227_GTACACCT-CATGAGGA_L001
KG228_ATGGCTGT-TTGACAGG_L001
KG231_GGAGGAAT-AACCGTTC_L001
KG234_GATCTTGC-TGCTCATG_L001
KG236_GCGTTAGA-GTATGCTG_L001
KG237_CAATAGCC-TGCGTAGA_L001
KG244_TTGCTGGA-TGTCCAAT_L001
KG246_GACTTGTG-TCGTGCAC_L001
KG249_CCTTCCAT-ATTGCGTG_L001
KG253_GTCATCGT-TTAAGCGG_L001
KG260_TCCATTGC-AGACCGTA_L001
KG263_TATGGCAC-TCTCGCAA_L001
KG264_AAGTGTGG-GATGGGCA_L001
KG265_AAGGAAGG-ACTCCGGT_L001
KG266_ACCAACAG-GAGCAAAC_L001
KG268_CCACAACA-TGTGGTAC_L001
KG270_GTCTGAGT-CACAAGTC_L001
KG271_ATTCGAGC-CAACCTAG_L001
KG274_TAGCACCT-GTTCCAGA_L001
KG275_ACGAATCC-TGGGTAAT_L001
KG276_GGAATGTC-GTGTTCCT_L001
KG278_GAACGGTT-ATTCCTCC_L001
KG281_TTGGGTAC-GTTGGAAG_L001
KG283_GACGTCAT-TCAACTGG_L001
KG288_CAGCTTCG-CAACGGGT_L001
KG290_CTCACAAC-GTCTCCTT_L001
KG291_CAACACAG-TGGCTATC_L001
KG296_GTACCACA-TGTTGTGG_L001
KG298_TGAGGACT-GTCTAGTG_L001
KG306_CTGAGCTC-GCTCTGAT_L001
KG307_AACCAGAG-AACCGAAG_L001
KG308_CAACTCCA-ACGGATTC_L001
KG313_AGCCTATC-GTTACGCA_L001
KG316_TACCGGCT-GCCAGCTA_L001
KG317_TCTTGTTT-GTGAAAGG_L001
KG323_AGGTAGGA-GGCGAACA_L001
KG325_AAGGAGAC-GTTGTGAG_L001
KG331_AGTCAGGT-GTGGTTAC_L001
KG333_GGACATCA-TGCAGGTA_L001
KG334_CATTATGG-ATGTGACT_L001
KG338_ACCTCTTC-TATGCAAG_L001
KG341_TCGCGCAA-TGCCTTGT_L001
KG343_CTCGAAAT-CTTACCTG_L001
KG344_TGTCACAC-CGAGATTA_L001
KG345_GCCTTAAC-AGCTCCTA_L001
KG346_CAACCGAG-CACAGGTA_L001
KG347_CAGGTAAG-ATTTCGAG_L001
KG348_AGGAACAC-GACATTCC_L001
KG349_CTTCGGTT-CTCTGGTT_L001
KG350_ACCGTAAG-TGGATATG_L001
KG351_TACGGTCT-GCAATGGA_L001
KG352_GGTTGAAC-AACGTGGA_L001
KG354_AAGACACC-TCGGTTAC_L001
KG356_ACGATATG-GAATTGTC_L001
KG360_GTCCTTGA-TCAGACGA_L001
KG365_CGTATCTC-GTACCTTG_L001
KG366_TCATCTCC-CTAGCAAG_L001
KG375_CAATGCGA-AGTTCGTC_L001
KG376_CCTAAACT-CTCCTCGT_L001
KG377_CGAGAGAA-TCTCTTCC_L001
KG378_CCTGTCAA-ACAGCCAT_L001
KG379_CCAGTTGA-ATGACGTC_L001
KG380_GTAACCAC-ACCTGACT_L001
KG381_TGATAGGC-ACTCAGAC_L001
KG382_AGTGACCT-CTCCTAGA_L001
KG383_TGGAAGCA-GTACTCTC_L001
KG384_CATATCCA-CTTACGGT_L001
KG385_GCAATTCC-TGTTCGAG_L001
KG386_GAATCACC-CTTCGTTC_L001
KG387_AAGCGACT-GGATTCGT_L001
KG388_CGCTAATC-ACGTCATA_L001
KG389_CTTCCAAC-GTACCCAA_L001
KG390_TTACGTGC-CAAGGTCT_L001
KG391_CCAGTATC-ATCTCGCT_L001
KG395_CACTGTAG-AGTCGCTT_L001
KG401_GGTATAGG-CTGGAGTA_L001
KG402_GTCAGTAT-AGGCTTCT_L001
KG404_CGCTCTTA-ACGGTATG_L001
KG406_GCACACAA-TCGTCAAG_L001
KG422_CACCTGTA-CTGACACA_L001
KG426_TTGCGAGA-GTGCCATA_L001
KG428_ATCCGCAG-AGTCGTGT_L001
KG430_AGCTAAGC-GTAACGAC_L001
KG433_GGCTTACT-AGGGAAAG_L001
KG434_CTTCACTG-GTCCTAAG_L001
KG435_CACGCAAT-ATGGAAGG_L001
KG438_TACTCCAG-CCTATACC_L001
KG441_CCAAGGTT-GCTATCCT_L001
KG447_CCGACTCT-CTGGTTCT_L001
KG454_CCTCGGGT-GTCGTGGC_L001
KG457_AGCTACCA-CTTAGGAC_L001
KG458_TACCTGCA-TGATGTCC_L001
KG460_CGAATACG-GTCGGTAA_L001
KG461_TCCTCATG-AGGTGTAC_L001
KG462_CTCAGAAG-AACTTGCC_L001
KG478_CGGCATTA-TGACTGAC_L001
KG480_CTTAGGAC-CAGTGAAG_L001
KG482_GAGGACCA-ACCGGCAT_L001
KG484_AGACCTTG-GCACGTAA_L001
KG486_AGACGCTA-TGTACACC_L001
KG488_GATATGAA-CTGCCGTA_L001
KG490_CATACCGT-TAAGAGCG_L001
KG492_CAGTGCTT-ACCTGGAA_L001
KG494_GTAACCGA-GGTGTCTT_L001
KG495_AAGGCGTA-GCGGATGG_L001
KG501_ATCAGAGC-GAGCTCAG_L001
KG502_GTCGTTAC-GCTTAGCT_L001
KG503_TCCCACGA-CTATTGAA_L001
KG508_GAGAGTAC-TGCTTCCA_L001
KG512_CAAGGTAC-GAGATACG_L001
KG513_ACGAGGAG-AGTTTAGG_L001
KG514_CTCGAACA-GGAATTGC_L001
KG515_CGAATTGC-CACCTTAC_L001
KG516_GTGCACGA-GCCTATCA_L001
KG517_CTTACAGC-CTTCACCA_L001
KG518_CAGCATAC-TCTAACGC_L001
KG523_AGTCACAT-CCATAATG_L001
KG524_CTTGCATA-GAAGAGGT_L001
KG526_GGCAAGTT-CTTCTGAG_L001
KG528_GACAATTC-CATATCGT_L001
KG529_ACACGACT-CTGCGGAT_L001
KG530_CTCGACTT-GGTCAGAT_L001
KG531_ATGTTCCT-GGCGCTGA_L001
KG533_TGTTCGCC-TCCTACCT_L001
KG535_ACAAGGCA-TTGCGCGA_L001
KG537_CCTTTCAC-AAACAAGA_L001
KG540_GGAAGAGA-TTCTCTCG_L001
KG542_TCAGCGCC-AGGAACAT_L001
KG544_CAGGTTCA-GGCGTTAT_L001
KG545_GCCAATCC-GATCGGAC_L001
KG546_AAGACCGT-CAATCGAC_L001
KG548_TCCACGTT-GTTCAACC_L001
KG552_ATGCCGGT-TGGTCCTC_L001
KG555_AGGATAGC-AACCTTGG_L001
KG557_GTCCTAAG-TGGTAGCT_L001
MSB285669_ACCCGTTG-CGAAGCTG_L001
MSB285671_TCTGGAAC-AGGTGCTA_L001
MSB285674_ACCGGAGT-CCTTCCTT_L001
MSB285675_GCCACGAC-ACCCGAGG_L001
MSB285676_AGAACCAG-AGAGTCGG_L001
MSB285680_ATTGGACA-TCCAGCAA_L001
MSB285682_GTTTGCTC-CTGTTGGT_L001
MSB285684_TTGCGTTA-TTCGGGAA_L001
MSB285685_TTCAATAG-TCGTGGGA_L001
MSB288874_CACTAGAC-AGTCCTCA_L001
MSB288917_CTAGGTTG-GCTCGAAT_L001
MSB289689_TGTGTCAG-TACAGGTG_L001
MSB293477_AGAAGCCT-ATACTGAC_L001
MSB293654_CTTGACGA-TTGTGTGC_L001
extra_barcode_GTCGATTG-ACGGTCTT_L001"

for i in $FILES
do
cutadapt -a "G{50}" -A "G{50}" -m 36 --cores=8 -o data/fastq/trim_polyG/"$i"_R1.polyGtrim.fastq.gz -p data/fastq/trim_polyG/"$i"_R2.polyGtrim.fastq.gz data/fastq/trimmed/"$i"_R1_001.trim.fastq.gz data/fastq/trimmed/"$i"_R2_001.trim.fastq.gz
done