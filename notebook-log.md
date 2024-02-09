# Data Description
These data are from my master's research on pika (*Ochotona*) evolutionary biology. I performed custom target enrichment on 412 nuclear candidate genes related to adaptive phenotypes. I also included unenriched libraries for low pass genome skimming to obtain mtDNA. Illumina sequencing (2x150) was performed on Novaseq6000, so my raw data are in fastq format.

# Data Processing Pipeline
## QC and Trimming
Demultiplexing and intial qc (with fastqc and multiqc) was performed by folks at the sequencing core.

I first trimmed adapters and low quality reads with trimmomatic (using default paramters outlined in trimmomatic vignette), then checked quality with fastqc and multiqc.

**Trimmomatic**
```
for infile in *_L001_R1_001.fastq.gz \
do \
base=$(basename ${infile} _L001_R1_001.fastq.gz) \
trimmomatic PE -threads 8 \
${base}_L001_R1_001.fastq.gz ${base}_L001_R2_001.fastq.gz \
trimmed/${base}_L001_R1_001.trim.fastq.gz untrim/${base}_L001_R1_001.untrim.fastq.gz \
trimmed/${base}_L001_R2_001.trim.fastq.gz untrim/${base}_L001_R2_001.untrim.fastq.gz \
ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE-2.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36 \
done
```
**fastqc**
```
for infile in *_001.trim.fastq.gz \
do \
base=$(basename ${infile} _001.trim.fastq.gz) \
fastqc -o fastqc --noextract ${base}_001.trim.fastq.gz \
done
```
**multiqc**
```
multiqc -o fastqc/multiqc fastqc
```

After examining fastqc files, I saw that some of my samples had poly-G tails. I removed these with cutadapt (see `scripts/polyGtrim.sh`) and checked with fastqc again.

**fastqc (after poly-G trim)**
```
for infile in *_001.polyGtrim.fastq.gz \
do \
base=$(basename ${infile} _001.polyGtrim.fastq.gz) \
fastqc -o trim_polyG/fastqc --noextract ${base}_001.polyGtrim.fastq.gz \
done
```
**multiqc (after poly-G trim)**
```
multiqc -o fastqc/multiqc fastqc
```

