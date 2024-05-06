# Data Description
These data are from my master's research on pika (*Ochotona*) evolutionary biology. I performed custom target enrichment on 412 nuclear candidate genes related to adaptive phenotypes. I also included unenriched libraries for low pass genome skimming to obtain mtDNA. Illumina sequencing (2x150) was performed on Novaseq6000, so my raw data are in fastq format. I sequenced 179 individuals, but am only using one representative for each non-*O. princeps* species (5 species) and one representaive for 37 *O. princeps* populations in my final analysis. 

# Data Processing Pipeline
## QC and Trimming
Demultiplexing and intial qc (with fastqc and multiqc) was performed by folks at the sequencing core.

I first trimmed adapters and low quality reads with trimmomatic (v0.39; using default paramters outlined in trimmomatic vignette), then checked quality with fastqc and multiqc.

**Trimmomatic**
```
for infile in *_L001_R1_001.fastq.gz 
do
base=$(basename ${infile} _L001_R1_001.fastq.gz)
java -jar /usr/share/java/trimmomatic.jar PE -threads 8 \
${base}_L001_R1_001.fastq.gz ${base}_L001_R2_001.fastq.gz \
trimmed/${base}_L001_R1_001.trim.fastq.gz untrim/${base}_L001_R1_001.untrim.fastq.gz \
trimmed/${base}_L001_R2_001.trim.fastq.gz untrim/${base}_L001_R2_001.untrim.fastq.gz \
ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE-2.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36 |& tee -a trimmomatic_output.txt
done 
```
**fastqc**
```
for infile in *_001.trim.fastq.gz 
do
base=$(basename ${infile} _001.trim.fastq.gz)
fastqc -o fastqc --noextract ${base}_001.trim.fastq.gz
done
```
**multiqc**
```
multiqc -o fastqc/multiqc fastqc
```

After examining fastqc files, I saw that some of my samples had poly-G tails. I removed these with cutadapt (v.1.18; see `scripts/polyGtrim.sh`) and checked with fastqc again.

**fastqc (after poly-G trim)**
```
for infile in *.polyGtrim.fastq.gz 
do 
base=$(basename ${infile} .polyGtrim.fastq.gz) 
fastqc -o fastqc --noextract ${base}.polyGtrim.fastq.gz
done
```
**multiqc (after poly-G trim)**
```
multiqc -o fastqc/multiqc fastqc
```

## Mapping
I first start by indexing the reference genome I used for creating my target capture baits (OchPri4.0) with bwa (v0.7.17).

```
bwa index GCF_014633375.1_OchPri4.0_genomic.fna.gz
```
Next I map my trimmed reads to the reference genome.
```
for infile in fastq/trim_polyG/*_L001_R1.polyGtrim.fastq.gz
do
base=$(basename ${infile} _L001_R1.polyGtrim.fastq.gz)
bwa mem -t 8 -M refgenome/GCF_014633375.1_OchPri4.0_genomic.fna.gz \
fastq/trim_polyG/${base}_L001_R1.polyGtrim.fastq.gz fastq/trim_polyG/${base}_L001_R2.polyGtrim.fastq.gz > sam/${base}.sam
done
```
I then convert sam to bam, sort bam files, and mark duplicates with samtools (v1.15.1).
```
for infile in *.sam
do
base=$(basename ${infile} .sam)
samtools fixmate -@ 8 -m -O bam ${base}.sam bam/${base}.bam
samtools sort -@ 8 -O bam -o bam/${base}.sorted.bam bam/${base}.bam
samtools markdup bam/${base}.sorted.bam bam/${base}.markdup.bam
done
```
## SNP Calling and Filtering
I then call SNPs with bcftools (v1.15.1), but must do this separately for mitochdonrial and nuclear genes due to the ploidy. 

I then reduced my sample list down to 1 high quality individual per lineage/species, drawing from different populations when possible for O. princeps

**Mitochondrial Genes**

Run the script `41_mitobam.sh`, which indexes bamfiles and extracts the mitochondrial coordinates from them, but only for the 41 samples that I want to use. 

Call variants and generate vcf. In mpileup I am disabling BAQ which helps to reduce false SNPs due to misalignments, applying minimum base quality of 20, and allowing for very high read depth as I will filter for this later. In call I am forcing a consensus for ploidy of 1 and only considering SNPs with p-value of 0.0001.
```
bcftools mpileup -B --min-BQ 20 -Ou --threads 8 -f refgenome/GCF_014633375.1_OchPri4.0_genomic.fna sam/bam/mitobam/*.markdup.mt.bam -d 8000 \
| bcftools call --threads 8 --ploidy 1 -p .0001 -mv -a GQ,GP -Oz -o 41_mtDNA.vcf.gz
```
Normalize vcf and rename samples.
```
bcftools norm -f refgenome/GCF_014633375.1_OchPri4.0_genomic.fna -Ov -o 41_mtDNA_norm.vcf 41_mtDNA.vcf.gz

bcftools reheader -s mt_vcf_rename_41.txt 41_mtDNA_norm.vcf -o 41_mtDNA_norm_renamed.vcf
```
Filter SNPs for minimum depth of 3 with vcftools (v0.1.16).
```
vcftools --vcf 41_mtDNA_norm_renamed.vcf --minDP 3 --recode --recode-INFO-all --out 41_mtDNA_norm_renamed_dp3
```
Then I filtered adjacent indels within 5bp (see this tutorial: https://samtools.github.io/bcftools/howtos/consensus-sequence.html)
```
bcftools filter --IndelGap 5 -Oz -o 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz 41_mtDNA_norm_renamed_dp3.recode.vcf
```
Index the vcf
```
tabix 41_mtDNA_norm_renamed_dp3_ig5.vcf.gz
```
___
**Nuclear Genes**

Generate a vcf with the `markdup.bam` files. I am using the same settings as the mtDNA vcf, except adding AD, DP, and SP annotations in case I want to filter for these settings later. I am also specifying a ploidy of 2. 

Move the 41 samples for the nuclear vcf into a new folder with `41_nucbam.sh` script.

```
bcftools mpileup -B --min-BQ 20 -a AD,DP,SP -Ou --threads 8 -f refgenome/GCF_014633375.1_OchPri4.0_genomic.fna sam/bam/41_samples/*.markdup.bam -d 2000 \
| bcftools call --threads 8 --ploidy 2 -p .0001 -mv -a GQ,GP -Oz -o 41.vcf.gz
```
___
Normalize vcf
```
bcftools norm -f refgenome/GCF_014633375.1_OchPri4.0_genomic.fna -Oz -o 41_norm.vcf.gz 41.vcf.gz
```

Filter SNPs for minimum depth of 4.
```
vcftools --gzvcf 41_norm.vcf.gz --minDP 4 --recode --recode-INFO-all --out 41_norm_dp4
```
Rename the samples
```
bcftools reheader -s vcf_rename.txt 41_norm_dp4.recode.vcf -o 41_norm_dp4_renamed.vcf
```
Then I filtered adjacent indels within 5bp (see this tutorial: https://samtools.github.io/bcftools/howtos/consensus-sequence.html)
```
bcftools filter --IndelGap 5 -Oz -o 41_norm_dp4_renamed_ig5.vcf.gz 41_norm_dp4_renamed.vcf
```
I am also trying a dataset where I remove indels completely, because I think most phylogenetic models don't use them anyways and they're causing me a lot of trouble with the consensus scripts. 
```
vcftools --gzvcf 41_norm_dp4_renamed_ig5.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | bgzip > 41_norm_dp4_renamed_ig5_noindel.vcf.gz
```
Index the vcf
```
tabix 41_norm_dp4_renamed_ig5.vcf.gz
tabix 41_norm_dp4_renamed_ig5_noindel.vcf.gz
```

## Generate fasta files for each gene
### Consensus Scripts
First I had to extract the cds coordinates for my nuclear candidate genes. I used `scripts/coordinates_minus_plus.sh` for this. Then I moved all of the coordiantes on the plus strand into a separate folder with grep, and manually moved the rest into another folder.

```
for file in *.txt; do
[[ $(grep "+" "$file")]] && mv "$file" plus
done
```

Use `scripts/mt_consensus.sh` to create fasta consensus files for each mitochondrial gene.  For plus strand nuclear genes, I run `scripts/nuc_consensus_plus.sh`. For minus strand genes, I run `scripts/nuc_consensus_minus.sh`. These scripts generate fasta's using the reference genome and variants in the vcf. I am using IUPAC codes for heterozygous sites, as RAxML can handle these, and coding missing data as "N". 

These scripts worked for most of the transcipts except for a couple dozen. I modified the scripts for these transcripts: `scripts/coordinates_product_sh` and `scripts/coordinates_mRNA.sh` because they didn't have the word "gene" in the sections I was extracting, but rather "mRNA" or "product". 

### Outgroup Homologous Genes
I use blast to get homologous genes from my outgroup (rabbit). I first make a blast database of the rabbit genome with BLAST+ (v2.12.0). 

```
makeblastdb -in GCF_000003625.3_OryCun2.0_genomic.fna -title OryCun2 -out OryCun2 -dbtype nucl -parse_seqids
```

Starting with nuclear genes, make fasta files for the pika genes
```
for infile in coordinates/all/*_coord.txt; do f=$(basename ${infile} _coord.txt); awk -F'\t' '{print $1}' coordinates/all/"$f"_coord.txt | samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna -r - > coordinates/gene_fasta/"$f"_ochpri4.fa; done
```
Concatenate these exons into a single-line fasta
```
for infile in coordinates/gene_fasta/*_ochpri4.fa; do f=$(basename ${infile} _ochpri4.fa); sed -s '1h;/>/d;H;$!d;x;s/\n//2g' coordinates/gene_fasta/"$f"_ochpri4.fa > coordinates/gene_fasta/"$f"_ochpri4_cat.fa; done
```
Run a blast search with these genes on the rabbit genome, keeping only the top subject hit
```
for infile in coordinates/gene_fasta/*_ochpri4_cat.fa; do f=$(basename ${infile} _ochpri4_cat.fa); blastn -query coordinates/gene_fasta/"$f"_ochpri4_cat.fa -db refgenome/Orycun2 -outfmt "6 sseq" -max_target_seqs 1 -task blastn -out OryCun2/"$f".fa 
done
```
Then concatenate each of the sequences from that subject to the same line, remove spaces, and add a fasta header
```
for infile in OryCun2/*.fa; do j=$(basename ${infile} .fa);
echo $(cat OryCun2/"$j".fa) | sed 's/ //g' | sed -e '1i\>Oryctolagus_cuniculus' > OryCun2/final_nuc/"$j".fa
done
```
Then add the rabbit sequence to the pikas fasta records
```
for infile in nuc_fasta/och_final/*_renamed.fa; do f=$(basename ${infile} _renamed.fa); cat nuc_fasta/och_final/"$f"_renamed.fa OryCun2/final_nuc/"$f".fa > nuc_fasta/och_ory/"$f".fa; done
```
___
Then mitochondrial genes
```
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:2745-3699 -o mt_ref_fasta/ND1.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:3908-4949 -o mt_ref_fasta/ND2.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:5333-6874 -o mt_ref_fasta/COX1.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:7018-7701 -o mt_ref_fasta/COX2.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:7773-7979 -o mt_ref_fasta/ATP8.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:7934-8613 -o mt_ref_fasta/ATP6.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:8614-9397 -o mt_ref_fasta/COX3.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:9468-9813 -o mt_ref_fasta/ND3.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:9883-10179 -o mt_ref_fasta/ND4L.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:10173-11550 -o mt_ref_fasta/ND4.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:11749-13560 -o mt_ref_fasta/ND5.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:13556-14080 -o mt_ref_fasta/ND6.fa
samtools faidx refgenome/GCF_014633375.1_OchPri4.0_genomic.fna NC_005358.1:14153-15290 -o mt_ref_fasta/CYTB.fa
```
Run a blast search with these genes on the rabbit genome, keeping only the top subject hit
```
for infile in mt_ref_fasta/*.fa; do f=$(basename ${infile} .fa); blastn -query mt_ref_fasta/"$f".fa -db refgenome/Orycun2 -outfmt "6 sseq" -max_target_seqs 1 -task blastn -out OryCun2/mt/"$f".fa 
done
```
Then concatenate each of the sequences from that subject to the same line, remove spaces, and add a fasta header
```
for infile in OryCun2/mt/*.fa; do j=$(basename ${infile} .fa);
echo $(cat OryCun2/mt/"$j".fa) | sed 's/ //g' | sed -e '1i\>Oryctolagus_cuniculus' > OryCun2/mt/final/"$j".fa
done
```
Then add the rabbit sequence to the pikas fasta records
```
for infile in mt_fasta_dp3/final/*_renamed.fa; do f=$(basename ${infile} _renamed.fa); cat mt_fasta_dp3/final/"$f"_renamed.fa OryCun2/mt/final/"$f".fa > mt_fasta_dp3/och_ory/"$f".fa; done
```

## Alignment
I then align with MAFFT (v7.490). MAFFT is unique in that it uses the Fast Fourier Transform for quickly finding homologous sequences, as it identifies correlations among physiochemical properties of amino acids. It can apply progressive alignment, iterative refinement, and structural alignment methods. The default is progressive alignment, which is what I use. MAFFT assumes that sequences are homologous with no genomic rearrangment, so it is limited by inversions, duplications, and translocations. 
I downloaded from this link: https://mafft.cbrc.jp/alignment/software/

For mtDNA
```
for infile in mt_fasta_dp3/och_ory/*.fa; do f=$(basename ${infile} .fa)
mafft mt_fasta_dp3/och_ory/"$f".fa > mt_fasta_dp3/och_ory/"$f".msa; done
```
For nDNA
```
for infile in nuc_fasta/och_ory/*.fa; do f=$(basename ${infile} .fa)
mafft nuc_fasta/och_ory/"$f".fa > nuc_fasta/och_ory/"$f".msa; done
```

However, the blast approach to obtaining my outgroup sequences made my alignment files super gappy, which made iqtree refuse to analyse some of the really bad ones. To attempt to get around this, I used MSA_trimmer to remove positions where more than half of the samples are missing data. https://github.com/LKremer/MSA_trimmer

MSA_trimmer is a bit old and requires an older version of biopython than what I had on my machine.
```
pip install biopython==1.77
```
```
for infile in mt_fasta_dp3/och_ory/*.msa; do f=$(basename ${infile} .msa)
python alignment_trimmer.py -c mt_fasta_dp3/och_ory/"$f".msa --trim_gappy 0.5; done
```
```
for infile in nuc_fasta/och_ory/*.msa; do f=$(basename ${infile} .msa)
python alignment_trimmer.py -c nuc_fasta/och_ory/"$f".msa --trim_gappy 0.5; done
```
I made a mistake while renaming the fastas, so rather than running consensus scripts all over again, I just changed the name here before building trees. But note that I fixed the rename files, so if I have to run consensus again, skip this step. 
```
for infile in nuc_fasta/och_ory/*_trimmed.msa; do f=$(basename ${infile} _trimmed.msa)
sed -i 's/SRM_UT_Mud/SN_UT_Mud/g' nuc_fasta/och_ory/"$f"_trimmed.msa; done
```
```
for infile in mt_fasta_dp3/och_ory/*_trimmed.msa; do f=$(basename ${infile} _trimmed.msa)
sed -i 's/SRM_UT_Mud/SN_UT_Mud/g' mt_fasta_dp3/och_ory/"$f"_trimmed.msa; done
```

# Tree estimation
## Maximum Likelihood
I first ran RAxML on my CTYB alignment (see below), it told me that I should change the site model from what I was using. Since I am building gene trees for 400+ genes and don't want to have to manually change my model for each gene, I am going to run this in IQ-tree (v2.0.7). IQ-tree is a maximum likelihood inference method that makes use of hill-climbing algorithms and stochastic purturbation. I chose IQ-tree because of it’s ModelFinder tool, fast bootstrap analysis, and positive feedback on accuracy. One review of IQ-tree found that 18% of runs were irreproducible, which may owe to their hill-climbing heuristic appraoch, however, I aim to address this by reporting my starting seed. IQ-tree assumes treelikeness among all of the sites, stationarity of site frequencies through time, homogeneity of substiution rates over time, and equal likelihood directionality among substituions. While IQ-tree allows for many user-selected options, I just set my bootstrap replicates to 1000 and my seed to 123. 

I installed IQtree from apt, as per these instructions: http://www.iqtree.org/doc/Quickstart
```
sudo apt-get install iqtree
```
For mtDNA genes, with 1000 bootstrap replicates, setting rabbit as my outgroup, allowing it to automatically set the number of threads, and a seed of 123. 
```
for infile in mt_fasta_dp3/och_ory/*_trimmed.msa; do f=$(basename ${infile} _trimmed.msa)
iqtree2 -s mt_fasta_dp3/och_ory/"$f"_trimmed.msa -bb 1000 -nt AUTO -T AUTO -seed 123 -o Oryctolagus_cuniculus --prefix "$f"_och_ory; done
```
For nDNA genes, with 1000 bootstrap replicates
```
for infile in nuc_fasta/och_ory/*_trimmed.msa; do f=$(basename ${infile} _trimmed.msa)
iqtree2 -s nuc_fasta/och_ory/"$f"_trimmed.msa -bb 1000 -nt AUTO -T AUTO -seed 123 -o Oryctolagus_cuniculus --prefix "$f"_och_ory; done
```
Then concatenate the trees into one file
```
cat *treefile > 429.tre
```

## Multispecies Coalescent
I am using Astral-III (v5.7.8) for my musltispecies coalescent approach. Astral is phylogenetic tree inference software that utilizes a multi-species coalescent appraoch for estimating an unrooted species tree from a list of unrooted gene trees. One strength is that Astral is very fast and by running a ML software on my genes and then Astral, there is a significant reduction in time than if I used a co-estimation approach. One limitation/weakness is that Astral does not generate a time-calibrated phylogeny. As far as assumptions, the final tree generated by Astral is only going to be as good as the input trees. Since I am using candidate genes that could be under selection, I will not be saying that this tree is the end-all species tree for my taxa. Since my gene trees are already generated, I will be running a pretty straightforward command. I will just tell it the input data, where to put the output, and have it save a log file. 

I installed the software following the guide here: https://github.com/smirarab/ASTRAL

```
java -jar Astral/astral.5.7.8.jar -i iqtrees/429.tre -o 429_astral.tre 2> 429_IQtree_astral.log
```

# Visualization and outlier detection
## Outlier detection
I conducted gene outlier detection with PhylteR (v0.9.11) in R (v4.4.0), which uses a multidimensional scaling approach to compare sets of three distance matrices. It removes outlier sequences and identifies genes that are still outliers after these edits are made. I input all 429 ML trees generated by IQ-tree with k=2. 

First I read in my ML gene trees with ape (v5.8).
```
install.packages(phylteR)
install.packages(ape)
install.packages(ggtree)
library(phylteR)
library(ape)
library(ggtree)
trees_429<-ape::read.tree("iqtrees/429.tre")
gene_list<-read.delim("gene_list.txt",header=FALSE)
```
Then I ran PhylteR and processed read out the results
```
phylter_results<-phylter(trees_429, gene.names=gene_list$V1)
phylter_results$Final$Outliers
write.phylter(phylter_results,"phylter.out")
```
## Tree visualization 
Code for visualizing the full 429 gene MSC tree with ggtree (v3.11.2)
```
full_tip_cat<-read.csv("tip_list.csv")
full_tip_cat<-as.data.frame(full_tip_cat)

astral_429<-read.tree("429_astral.tre")
rooted_429<-ape::root(astral_429,outgroup = "Oryctolagus_cuniculus",resolve.root=TRUE)

astral_429_tree<-ggtree(rooted_429)  +
  geom_text2(aes(label=label, subset = !is.na(as.numeric(label))),nudge_x=-.23,nudge_y=.35) +
  ggplot2::xlim(0, 22)

astral_429_tree %<+% full_tip_cat + 
  geom_tiplab(aes(fill=factor(Cat)),
              geom="label",size=3.2,alpha=0.7) + 
  scale_fill_manual("Taxa Category",values=c("SRM"="#FF0000",
                                             "CU"="#FFA500",
                                             "CR"="#33a02c",
                                             "SN"="#990099",
                                             "NRM_NW"="#0066CC",
                                             "NRM_SE"="#999999",
                                             "Pika_Asia"="mediumturquoise",
                                             "Pika_NorthAmerica"="lightsteelblue2",
                                             "Ochotona_Asia"="yellow2",
                                             "Outgroup"="indianred")) +
  theme(legend.position=c(0.85,0.7)) + geom_treescale()
```
Code for the 416 nuclear gene tree
```
astral_416<-read.tree("416_astra.tre")
rooted_416<-ape::root(astral_416,outgroup = "Oryctolagus_cuniculus",resolve.root = TRUE)

astral_416_tree<-ggtree(rooted_416)  +
  geom_text2(aes(label=label, subset = !is.na(as.numeric(label))),nudge_x=-.33,nudge_y=.35) +
  ggplot2::xlim(0, 18)

astral_nuc_tangle<-astral_416_tree %<+% full_tip_cat + 
  geom_tiplab(aes(fill=factor(Cat)),
              geom="label",size=3.2,alpha=0.7) + 
  scale_fill_manual("Taxa Category",values=c("SRM"="#FF0000",
                                             "CU"="#FFA500",
                                             "CR"="#33a02c",
                                             "SN"="#990099",
                                             "NRM_NW"="#0066CC",
                                             "NRM_SE"="#999999",
                                                       "Pika_Asia"="mediumturquoise",
                                                       "Pika_NorthAmerica"="lightsteelblue2",
                                                       "Ochotona_Asia"="yellow2",
                                                       "Outgroup"="indianred")) +
  theme(legend.position = "none") + geom_treescale()
```
Code for the mtDNA tree
```
astral_mt<-read.tree("mt_astral.tre")
rooted_mt<-ape::root(astral_mt,outgroup="Oryctolagus_cuniculus")

astral_mt_tree<-ggtree(rooted_mt) + geom_text2(aes(label=label, subset = !is.na(as.numeric(label))),nudge_x=0.4,nudge_y=.32)

astral_mt_tangle<-astral_mt_tree %<+% full_tip_cat + 
  geom_tiplab(aes(fill=factor(Cat)), hjust=1,
              geom="label",size=3.2,alpha=0.7) + 
  scale_fill_manual("Taxa Category",values=c("SRM"="#FF0000",
                                             "CU"="#FFA500",
                                             "CR"="#33a02c",
                                             "SN"="#990099",
                                             "NRM_NW"="#0066CC",
                                             "NRM_SE"="#999999",
                                             "Pika_Asia"="mediumturquoise",
                                             "Pika_NorthAmerica"="lightsteelblue2",
                                             "Ochotona_Asia"="yellow2",
                                             "Outgroup"="indianred")) +
  theme(legend.position = "none") + scale_x_reverse(limit=c(18,0))  + geom_treescale(x=-9)
```
Then to get the nDNA and mtDNA trees side-by-side. Eventaully I would like to make a tanglegram, but I was having trouble with the code for now. 
```
library(cowplot)
cowplot::plot_grid(astral_nuc_tangle, astral_mt_tangle, ncol=2)
```
Finally, code for the outlier gene trees
```
outlier_trees<-lapply(Sys.glob("iqtrees/429_outliers/*.treefile"),read.tree)

rooted_tree<-list()
for (i in (1:24)) {
  rooted_tree[[i]]<-root(outlier_trees[[i]],outgroup = "Oryctolagus_cuniculus",resolve.root=TRUE)
}

outlier_tree<-list()
for (i in (1:24)) {
  outlier_tree[[i]]<-ggtree(rooted_tree[[i]]) + 
    geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70)) +
    ggplot2::xlim(0, 0.35) 
}

color_tree<-list()
for (i in (1:24)) {
  color_tree[[i]]<-outlier_tree[[i]] %<+% full_tip_cat + 
    geom_tiplab(aes(fill=factor(Cat)),
                geom="label", alpha=0.7,offset=0.004) + 
    scale_fill_manual("Subgenus_Region_Species",values=c("SRM"="#0066CC",
                                                         "CU"="#0066CC",
                                                         "CR"="#0066CC",
                                                         "SN"="#0066CC",
                                                         "NRM_NW"="#0066CC",
                                                         "NRM_SE"="#0066CC",
                                                         "Pika_Asia"="mediumturquoise",
                                                         "Pika_NorthAmerica"="lightsteelblue2",
                                                         "Ochotona_Asia"="yellow2",
                                                         "Outgroup"="indianred")) +
    theme(legend.position.inside=c(0.7,0.7)) + geom_treescale()

plot_list(color_tree[[1]],color_tree[[6]])
plot_list(color_tree[[10]],color_tree[[13]])
```
___
# Extra Items from the Course
These are extra items that were for the course homeworks, but I didn't end up using for my final project.

## Distance and parsimony
### Distance
For a distance-based tree, I used the ape package in R and used the default model of evolution (K80). The K80 model assumes that all substitutions have the same probability except that transitions and transversions are different from one another. I then used the distance matrix to generate a neighbor-joining tree, which is unrooted and does not assume a constant molecular clock. 

The advantage of using distance-based methods is that it is fast and scalable to large datasets, as it can generate a tree without having to search the entire tree space. The disavantage is that the tree is an approximation fo the optimum, as it has not searched the entire tree space.

In addition to setting the model of evolution, I can set a gamma parameter (defult is no gamma), set how to handle missing data (default is to delete sites with at least one sample with missing data), set to calculate variances from the distnaces (defualt is no), how the base frequencies should be calculatd (default is from the entire sequence), and whether to return as a matrix. I used all default parameters as I am only using this as a comparison to other methods. 

Below I am providing an example for generating a distance-based tree with my CYTB alignment. 
```
library(ape)
library(adegenet)

dna <- fasta2DNAbin(file=CYTB.msa)
D <- dist.dna(dna)
tre <- nj(D)
tre <- ladderize(tre)
plot(tre, cex=.6)
title("NJ Distance Tree for CYTB")
```
### Parsimony
For a parsimony tree, I am using the phanghorn package in R. I start by converting the alignment into phydat format and generate a distance-based neighbor-joining tree to start the NNI search with, using the "raw" model of evolution, which is just based on the number of sites that differ in each pair of sequences. I then generate the parsimony score and generate the parsimony tree with the phydat object and the starting distance tree.

The main assumption of parsimony is independence among sites. Parsimony is not model-based and is useful when there are computational limitations with model-based methods. The disadvantage is that it sometimes produces inconsistent trees and cannot handle pathological inequalities effectively.

In phanghorn, I can change the method from fitch (default) to sankoff and whether to use NNI (default) or SPR rearrangements. I again use the default settings because I am using it as a comparison for other methods. 

Below I am providing an example for generating a parsimony tree with my CYTB alignment.
```
library(adegenet)
library(phanghorn)

dna <- fasta2DNAbin(file=CYTB.msa)
dna2 <- as.phyDat(dna)
tre.ini <- nj(dist.dna(dna,model="raw"))
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2)
plot(tre.pars, cex=0.6)
title("Parsimony Tree for CYTB")
```
## Maximum Likelihood
I first tried RAxML to generate maximum likelihood trees. RAxML is self-described as a fast maximum likelihood tree search algorithm that returns trees with good likelihood scores. Strengths of RAxML are that it supports a wide variety of input data types, offers parallelization, and can automtically compute how many bootstrap values are necessary for each dataset. When compared with another leading ML software (IQTree), IQTree may require fewer replicates to find the best tree as it had lower variance than RAxML in certain test runs. I use the high performance computing (HPC) version and elected to use the GTRGAMMA model. The GTR model assumes that mutation process is the same at every branch of the tree, sites evolve independently, and all sites evolve the same rate. The last assumption is often violated in real data, so I allow substitution rates to vary at each site according to a gamma distrubtion. 

You can specify many options, but I use `-T` to set threads, `-m` to set the model, `-p` to set the seed for parsimony inferencs, `-#` to set the number of trees to compute, `-s` to set the input sequence, and `-n` to set the output name.

Below I am providing an example for generating a ML tree  with boostrap replicates using my CYTB alignment.
___

Start with computing 100 ML trees, which are all derived from distinct starting trees. The best supported tree will be called RAxML_bestTree.CYTB. I am using the GTRGAMMA model, running with 4 threads, and setting a seed of 12345.
```
raxmlHPC-PTHREADS -T 4 -m GTRGAMMA -p 12345 -# 100 -s CYTB.msa -n CYTB
```
To get bootstrap values I include the `-b` flag, also using the 12345 seed. I allow RAxML to choose the number of bootstrap replicates needed with `-# autoMRE`. This stands for extended majority-rule consensus tree criterion. 
```
raxmlHPC-PTHREADS -T 4 -m GTRGAMMA -p 12345 -b 12345 -# autoMRE  -s CYTB.msa -n CYTB_bs
```
Finally, I use the bootstrap file to draw bipartitions on the best ML tree. 
```
raxmlHPC-PTHREADS -T 4 -m GTRGAMMA -p 12345 -f b -t RAxML_bestTree.CYTB -z RAxML_bootstrap.CYTB_bs -n CYTB_bp
```

## Bayesian 
I am using MrBayes to generate a bayesian comparison for my CYTB alignment. MyBayes is a commonly used Bayesian Inference software that uses a variant of Markov chain Monte Carolo. Some major strengths of MrBayes are that it is easy to install and run and is widely used, so looking for information on paramters for various datasets can be easily acheived. A weakness of MrBayes, and most Bayesian analyses, is that there is rarely a priori for information for setting priors on trees and parameters. The independent exponential prior on branch length has been known to produce overly long branch trees with some datasets (Yang and Rannala 2012). MrBayes also assumes that aligned sequences are without errors and lacks ability to add new evolutionary models at the speed that RevBayes is able to. User choices include the number of generations, how often to save parameter values and tree topologies, how many independent runs, the number of chains, evolutionary model parameters, and priors for branch lengths, gamma shape paramter, kappa, and base frequencies. Assumptions will depend on the model of evolution, but in general, it is assumed that the mutation process is the same at every branch of the tree, sites evolve independently, and all sites evolve the same rate. The last assumption is often violated in real data, so I allow substitution rates to vary at each site according to a gamma distrubtion. 
___
I first start by converting my fasta to nexus using seqmagick. Also included is the code for how to install seqmagick.
```
pip install seqmagick
seqmagick convert --output-format nexus --alphabet dna CYTB.msa CYTB.nex
```
I then add this block to the nexus file. I am running with 10,000,000 generations, set the sample frequency to 1000, and number of chains to 5. Here I am using the same `lset` settings that were used in the tutorial (2-parameter substitution matrix and to have rates vary across sites along a gamma distrubtion with 4 categories). I am also using the same priors from the tutorial for now, as I have been having trouble finding priors specified for other Ochotona studies that used MrBayes in the literature. This includes an unconstrained branch length prior with an exponential mean of 1/10, the shape paramter is exponential with a mean of 1.0, the kappa prior with a beta(1,1) distrubtion, and a flat Dirichlet distribution. 
```
begin mrbayes;
 set autoclose=yes;
 prset brlenspr=unconstrained:exp(10.0);
 prset shapepr=exp(1.0);
 prset tratiopr=beta(1.0,1.0);
 prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
 lset nst=2 rates=gamma ngammacat=4;
 mcmcp ngen=10000000 samplefreq=1000 printfreq=1000 nruns=1 nchains=5 savebrlens=yes;
 outgroup dauurica;
 mcmc;
 sumt;
end;
```
I installed MrBayes with homebrew and ran with the following commands
```
brew install mrbayes
mb CYTB.nex
```