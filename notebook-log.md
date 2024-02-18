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

## Generate fasta files for each gene and align with MAFFT
Use `scripts/mt_consensus.sh` to create fasta consensus files for each mitochondrial gene and `scripts/nuc_consensus.sh` for candidate nuclear genes. These scripts generate fasta's using the reference genome and variants in the vcf. I am using IUPAC codes for heterozygous sites, as RAxML can handle these, and coding missing data as "N". 

Before running these two scripts, I activate mafft
```
conda activate mafft
```

## Align
The gene fasta's are probably already aligned, but to be sure, I realigned with MAFFT. Note that this was already performed in `mt_consensus.sh` and `nuc_consensus.sh`, but this was the command:

```
mafft {gene_name}.fa > {gene_name}.msa
```
