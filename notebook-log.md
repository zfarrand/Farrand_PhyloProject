# Data Description
These data are from my master's research on pika (*Ochotona*) evolutionary biology. I performed custom target enrichment on 412 nuclear candidate genes related to adaptive phenotypes. I also included unenriched libraries for low pass genome skimming to obtain mtDNA. Illumina sequencing (2x150) was performed on Novaseq6000, so my raw data are in fastq format.

# Data Processing Pipeline
## QC and Trimming
Demultiplexing and intial qc (with fastqc and multiqc) was performed by folks at the sequencing core.

I first trimmed adapters and low quality reads with trimmomatic (using default paramters outlined in trimmomatic vignette), then checked quality with fastqc and multiqc.

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

After examining fastqc files, I saw that some of my samples had poly-G tails. I removed these with cutadapt (see `scripts/polyGtrim.sh`) and checked with fastqc again.

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
I first start by indexing the reference genome I used for creating my target capture baits (OchPri4.0) with bwa.

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
I then convert sam to bam, sort bam files, and mark duplicates.
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
I then call SNPs with bcftools, but must do this separately for mitochdonrial and nuclear genes due to the ploidy. 

I then reduced my sample list down to 3 high quality individuals per lineage/species, drawing from different populations when possible. 

**Mitochondrial Genes**

Start by indexing bam files.
```
for infile in *.markdup.bam
do
base=$(basename ${infile} .markdup.bam)
samtools index ${base}.markdup.bam
done
```
Then pull out the mtDNA sequences based on coordinates from the reference genome.
```
for infile in *.markdup.bam
do
base=$(basename ${infile} .markdup.bam)
samtools view -@ 8 -o mitobam/${base}.mt.bam ${base}.markdup.bam NC_005358.1:1-16481
done
```
Call variants and generate vcf. In mpileup I am disabling BAQ which helps to reduce false SNPs due to misalignments, skipping indels, applying minimum base quality of 20, and allowing for very high read depth as I will filter for this later. In call I am forcing a consensus for ploidy of 1 and only considering SNPs with p-value of 0.0001.
```
bcftools mpileup -B -I --min-BQ 20 -Ou --threads 8 -f GCF_014633375.1_OchPri4.0_genomic.fna *.markdup.bam -d 8000 \
| bcftools call --threads 8 --ploidy 1 -p .0001 -mv -a GQ,GP -Oz -o mtDNA.vcf.gz
```
Normalize vcf, limit samples to one high quality individual per location, and rename samples.
```
bcftools norm -f GCF_014633375.1_OchPri4.0_genomic.fna -Oz -o mtDNA_norm.vcf.gz mtDNA.vcf.gz

vcftools --gzvcf mtDNA_norm.vcf.gz --remove mt_remove.txt --recode --recode-INFO-all --out mt_46.vcf

bcftools reheader -s mt_46_rename.txt mt_33.recode.vcf -Oz -o mt_46_renamed.vcf.gz
```
Filter SNPs for minimum depth of 3.
```
vcftools --gzvcf mt_46_renamed.vcf.gz --minDP 3 --recode --recode-INFO-all --out mt_46_dp3
```
Compress and index the vcf
```
bgzip mt_46_dp3.recode.vcf
tabix mt_46_dp3.recode.vcf.gz
```
**I am still working on the nuclear vcf scripts as of 2-19**

## Generate fasta files for each gene and align with MAFFT
### Consensus Scripts
First I had to extract the cds coordinates for my nuclear candidate genes. I used `scripts/coordinates.sh` for this. Then I moved all of the coordiantes on the plus strand into a separate folder with grep, and manually moved the rest into another folder.

```
for file in *.txt; do
[[ $(grep "+" "$file")]] && mv "$file" plus
done
```


Use `scripts/mt_consensus.sh` to create fasta consensus files for each mitochondrial gene and `scripts/nuc_consensus.sh` for candidate nuclear genes. These scripts generate fasta's using the reference genome and variants in the vcf. I am using IUPAC codes for heterozygous sites, as RAxML can handle these, and coding missing data as "N". 

**Note that I am still working on the consensus scripts as of 2-19**

Before running these two scripts, I activate mafft
```
conda activate mafft
```

### Alignment
I then realign with MAFFT. MAFFT is unique in that it uses the Fast Fourier Transform for quickly finding homologous sequences, as it identifies correlations among physiochemical properties of amino acids. It can apply progressive alignment, iterative refinement, and structural alignment methods. The default is progressive alignment, which is what I use. MAFFT assumes that sequences are homologous with no genomic rearrangment, so it is limited by inversions, duplications, and translocations. 

Note that this was already performed in `mt_consensus.sh` and `nuc_consensus.sh`, but this was the command:

```
mafft {gene_name}.fa > {gene_name}.msa
```
# Tree building
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
I am using RAxML to generate maximum likelihood trees. RAxML is self-described as a fast maximum likelihood tree search algorithm that returns trees with good likelihood scores. Strengths of RAxML are that it supports a wide variety of input data types, offers parallelization, and can automtically compute how many bootstrap values are necessary for each dataset. When compared with another leading ML software (IQTree), IQTree may require fewer replicates to find the best tree as it had lower variance than RAxML in certain test runs. I use the high performance computing (HPC) version and elected to use the GTRGAMMA model. The GTR model assumes that mutation process is the same at every branch of the tree, sites evolve independently, and all sites evolve the same rate. The last assumption is often violated in real data, so I allow substitution rates to vary at each site according to a gamma distrubtion. 

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