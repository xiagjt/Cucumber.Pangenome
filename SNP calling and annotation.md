# SNP Calling and Annotation Pipeline using NGS Data

## Overview
Comprehensive pipeline for calling single nucleotide polymorphisms (SNPs) and small insertions/deletions (InDels) from next-generation sequencing (NGS) data using GATK best practices workflow, followed by functional annotation using ANNOVAR.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Step 0: Software Installation](#step-0-software-installation)
- [Step 1: Read Mapping](#step-1-read-mapping)
- [Step 2: Variant Calling (Per Sample)](#step-2-variant-calling-per-sample)
- [Step 3: Joint Genotyping](#step-3-joint-genotyping)
- [Step 4: Variant Filtering](#step-4-variant-filtering)
- [Step 5: Functional Annotation](#step-5-functional-annotation)


---

## Prerequisites

### Required Data
- **Reference genome**: CLv4.0 (CLv4.fa)
- **Gene annotation**: CLv4.gtf
- **NGS reads**: Paired-end FASTQ files (cleaned/filtered)

### Required Software
- BWA (Burrows-Wheeler Aligner)
- SAMtools
- GATK (Genome Analysis Toolkit)
- ANNOVAR
- VCFtools
- Perl (for ANNOVAR scripts)

---

## Step 0: Software Installation

### Install via Conda (Recommended)

```bash
# Create a new conda environment
conda create -n snp_calling python=3.8

# Activate environment
conda activate snp_calling

# Install required software
conda install -c bioconda bwa samtools gatk4 vcftools

# ANNOVAR needs to be installed separately (requires registration)
# Download from: http://www.openbioinformatics.org/annovar/
```

### Verify Installation

```bash
# Check BWA
bwa

# Check SAMtools
samtools --version

# Check GATK
gatk --version

# Check VCFtools
vcftools --version

# Check ANNOVAR
perl table_annovar.pl --help
```

---

## Step 1: Read Mapping

### Overview
Map clean NGS reads to reference genome using BWA-MEM algorithm, then sort and index BAM files.

### 1.1 Index Reference Genome

```bash
# Create BWA index (only need to do once)
bwa index CLv4.fa
```

**Output Files:**
- `CLv4.fa.amb`
- `CLv4.fa.ann`
- `CLv4.fa.bwt`
- `CLv4.fa.pac`
- `CLv4.fa.sa`

**Time:** ~30 minutes for 350 Mb genome

### 1.2 Map Reads for Each Sample

```bash
# Map reads with read group information
bwa mem \
  -t 20 \
  -R '@RG\tID:RXXXID\tSM:XXX' \
  CLv4.fa \
  sample_R1.fq \
  sample_R2.fq | \
  samtools view -@ 20 -bh - | \
  samtools sort -@ 20 -T sample_tmp -o sample.bam -

# Index BAM file
samtools index sample.bam
```

### Parameters Explained

**BWA mem:**
- `-t 20`: Number of threads
- `-R`: Read group information (REQUIRED for GATK)
  - `ID`: Read group identifier (unique ID, can be sample name or lane ID)
  - `SM`: Sample name (IMPORTANT: use consistent sample name across all analyses)

**SAMtools view:**
- `-@ 20`: Number of threads
- `-b`: Output BAM format
- `-h`: Include header

**SAMtools sort:**
- `-@ 20`: Number of threads
- `-T sample_tmp`: Temporary file prefix
- `-o`: Output file name

### Read Group Information

**Critical Importance:**
Read group (@RG) tags are REQUIRED for GATK. The minimum required tags are:
- `ID`: Read group identifier
- `SM`: Sample name (must be consistent for multi-sample calling)

**Example Read Groups:**

```bash
# Minimal read group (as shown in methods)
-R '@RG\tID:Sample001\tSM:Sample001'

# Full read group (optional but recommended)
-R '@RG\tID:Sample001_L1\tSM:Sample001\tPL:ILLUMINA\tLB:Lib1'
```

### Quality Checks

```bash
# Check mapping statistics
samtools flagstat sample.bam

# Check coverage
samtools depth sample.bam | \
  awk '{sum+=$3} END {print "Average depth:", sum/NR}'

# Check read groups
samtools view -H sample.bam | grep '^@RG'
```

**Expected Results:**
- Mapping rate: >95% for good quality data
- Average depth: 10-30× for resequencing
- Read groups: Should match input -R parameter

### Batch Processing Script

```bash
#!/bin/bash
# Script: map_all_samples.sh

ref="CLv4.fa"
threads=20

# Loop through all sample pairs
for fq1 in *_R1.fq; do
  fq2=${fq1/_R1.fq/_R2.fq}
  sample=$(basename $fq1 _R1.fq)
  
  echo "Processing: $sample"
  
  bwa mem \
    -t $threads \
    -R "@RG\tID:${sample}\tSM:${sample}" \
    $ref $fq1 $fq2 | \
    samtools view -@ $threads -bh - | \
    samtools sort -@ $threads -T ${sample}_tmp -o ${sample}.bam -
  
  samtools index ${sample}.bam
  
  echo "Completed: $sample"
done
```

---

## Step 2: Variant Calling (Per Sample)

### Overview
Call variants for each sample individually using GATK HaplotypeCaller in GVCF mode. This enables efficient joint genotyping across many samples.

### 2.1 Call Variants per Sample

```bash
gatk HaplotypeCaller \
  -R CLv4.fa \
  --emit-ref-confidence GVCF \
  -I sample.bam \
  -O sample.gvcf
```

### Parameters Explained

- `-R CLv4.fa`: Reference genome
- `--emit-ref-confidence GVCF`: Output GVCF format (genomic VCF)
  - Records both variant and reference confidence
  - Essential for joint genotyping
- `-I`: Input BAM file
- `-O`: Output GVCF file

### GVCF Format

**What is GVCF?**
- Genomic VCF: Contains both variant and reference blocks
- Stores genotype likelihoods for all sites
- Enables efficient joint calling across large cohorts
- Smaller than storing raw BAM files

**Example GVCF Records:**
```
##fileformat=VCFv4.2
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample001
Chr1    100     .       A       .       .       .       END=200 GT:DP:GQ:MIN_DP:PL      0/0:25:72:20:0,60,900
Chr1    201     .       G       A       150.5   .       .       GT:AD:DP:GQ:PL          0/1:10,15:25:99:180,0,150
```

### Batch Processing Script

```bash
#!/bin/bash
# Script: call_variants_all.sh

ref="CLv4.fa"

for bam in *.bam; do
  sample=$(basename $bam .bam)
  echo "Calling variants for: $sample"
  
  gatk HaplotypeCaller \
    -R $ref \
    --emit-ref-confidence GVCF \
    -I $bam \
    -O ${sample}.gvcf
  
  echo "Completed: $sample"
done
```

### Expected Runtime
- **Single sample:** 2-6 hours (30× coverage, 350 Mb genome)
- **Depends on:**
  - Sequencing depth
  - Genome size
  - CPU/memory resources

---

## Step 3: Joint Genotyping

### Overview
Combine GVCFs from all samples and perform joint genotyping to call variants across the entire cohort.

### 3.1 Combine GVCFs

```bash
# Combine all GVCF files
gatk CombineGVCFs \
  -R CLv4.fa \
  -V sample1.gvcf \
  -V sample2.gvcf \
  -V sample3.gvcf \
  -O cohort.combined.gvcf
```

**Note:** For many samples, you can list all GVCF files:
```bash
# Create list of all GVCFs
ls *.gvcf > gvcf_list.txt

# Use -V multiple times or from file
for gvcf in $(cat gvcf_list.txt); do
  echo "-V $gvcf"
done
```

### 3.2 Joint Genotyping

```bash
gatk GenotypeGVCFs \
  -R CLv4.fa \
  -V cohort.combined.gvcf \
  -O cohort.vcf
```

### Parameters Explained

**CombineGVCFs:**
- `-R`: Reference genome
- `-V`: Input GVCF file(s) - specify multiple times
- `-O`: Output combined GVCF

**GenotypeGVCFs:**
- `-R`: Reference genome
- `-V`: Combined GVCF file
- `-O`: Output raw VCF with genotypes

### Output VCF Format

**Raw VCF contains:**
- All detected variants (SNPs and InDels)
- Genotypes for all samples
- Quality scores and annotations
- Both filtered and unfiltered variants

**Example VCF Record:**
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample001       Sample002       Sample003
Chr1    1234    .       A       G       250.5   .       AC=3;AF=0.5     GT:AD:DP:GQ:PL  0/1:10,12:22:99:180,0,150       1/1:0,25:25:75:450,75,0 0/0:30,0:30:90:0,90,1350
```

---

## Step 4: Variant Filtering

### Overview
Filter variants based on minor allele frequency (MAF) and missing data rate to obtain high-quality SNP set.

### 4.1 Filtering with VCFtools

```bash
vcftools \
  --vcf cohort.vcf \
  --recode \
  --maf 0.05 \
  --max-missing 0.5 \
  --out filter.vcf
```

### Parameters Explained

- `--vcf`: Input VCF file
- `--recode`: Write new VCF file
- `--maf 0.05`: Minimum minor allele frequency (5%)
  - Removes rare variants
  - Focus on common polymorphisms
- `--max-missing 0.5`: Maximum missing data rate (50%)
  - Keep sites with ≤50% missing genotypes
  - At least 50% of samples must have data
- `--out`: Output prefix

### Output Files

- `filter.vcf.recode.vcf` - Filtered VCF file
- `filter.vcf.log` - Filtering statistics

### Advanced Filtering Options

```bash
# More comprehensive filtering
vcftools \
  --vcf cohort.vcf \
  --recode \
  --recode-INFO-all \
  --maf 0.05 \
  --max-missing 0.5 \
  --minQ 30 \
  --min-meanDP 10 \
  --max-meanDP 100 \
  --minDP 5 \
  --remove-indels \
  --out filter.vcf.snps
```

**Additional Parameters:**
- `--recode-INFO-all`: Keep all INFO fields
- `--minQ 30`: Minimum quality score
- `--min-meanDP 10`: Minimum mean depth across samples
- `--max-meanDP 100`: Maximum mean depth (avoid repeat regions)
- `--minDP 5`: Minimum depth per sample
- `--remove-indels`: Keep only SNPs (exclude InDels)

### Filtering Statistics

```bash
# Count variants before filtering
grep -v '^#' cohort.vcf | wc -l

# Count variants after filtering
grep -v '^#' filter.vcf.recode.vcf | wc -l
```

### Typical Filtering Results

**Before Filtering:**
- Total variants: 2,000,000 - 5,000,000
- SNPs: 80-85%
- InDels: 15-20%

**After Filtering (MAF > 0.05, Missing < 0.5):**
- Retained variants: 500,000 - 1,500,000
- Retention rate: 25-50%
- High-quality common SNPs

---

## Step 5: Functional Annotation

### Overview
Annotate filtered SNPs with gene information, coding consequences, and functional predictions using ANNOVAR.

### 5.1 Prepare Reference Files

```bash
# Convert GTF to genePred format
gtfToGenePred -genePredExt CLv4.gtf CLv4_refGene.txt

# Extract transcript sequences
perl retrieve_seq_from_fasta.pl \
  --format refGene \
  --seqfile CLv4.fa \
  CLv4_refGene.txt \
  --out CLv4_refGeneMrna.fa
```

**Required Files:**
- `CLv4_refGene.txt` - Gene structure in genePred format
- `CLv4_refGeneMrna.fa` - Transcript sequences

**Note:** Only need to prepare these files once per reference genome.

### 5.2 Annotate Variants

```bash
table_annovar.pl \
  filter.vcf \
  ./ \
  -buildver CLv4 \
  -out prefix \
  -remove \
  -protocol refGene \
  -operation g \
  -nastring . \
  -vcfinput
```

### Parameters Explained

- `filter.vcf`: Input VCF file (filtered variants)
- `./`: Database directory (contains CLv4_refGene.txt and CLv4_refGeneMrna.fa)
- `-buildver CLv4`: Genome build version (must match file prefix)
- `-out prefix`: Output file prefix
- `-remove`: Remove temporary files after completion
- `-protocol refGene`: Annotation protocol (gene-based annotation)
- `-operation g`: Operation type (g = gene-based)
- `-nastring .`: String for missing values
- `-vcfinput`: Input is VCF format

### Output Files

**Main Outputs:**
- `prefix.CLv4_multianno.txt` - Tab-delimited annotation table
- `prefix.CLv4_multianno.vcf` - Annotated VCF file

**Annotation Table Format:**
```
Chr     Start   End     Ref     Alt     Func.refGene    Gene.refGene    GeneDetail.refGene      ExonicFunc.refGene      AAChange.refGene
Chr1    12345   12345   A       G       exonic          CsaV4_1G000010  .                       synonymous SNV          CsaV4_1G000010:c.123A>G:p.Ala41Ala
Chr1    23456   23456   C       T       exonic          CsaV4_1G000020  .                       nonsynonymous SNV       CsaV4_1G000020:c.456C>T:p.Pro152Leu
Chr1    34567   34567   G       A       intronic        CsaV4_1G000030  .                       .                       .
Chr1    45678   45678   T       C       upstream        CsaV4_1G000040  dist=2345               .                       .
```

**Annotation Fields:**
- **Func.refGene**: Functional location
  - exonic, intronic, UTR5, UTR3, splicing, upstream, downstream, intergenic
- **Gene.refGene**: Gene name/ID
- **ExonicFunc.refGene**: Exonic consequence
  - synonymous SNV, nonsynonymous SNV, stopgain, stoploss, frameshift insertion, frameshift deletion
- **AAChange.refGene**: Amino acid change
  - Format: Gene:c.DNA_change:p.Protein_change

### Summarize Annotation Results

```bash
# Count variants by functional category
cut -f6 prefix.CLv4_multianno.txt | tail -n +2 | sort | uniq -c

# Count exonic variants by type
awk 'NR>1 && $6=="exonic" {print $9}' prefix.CLv4_multianno.txt | sort | uniq -c

# Extract nonsynonymous SNVs
awk 'NR>1 && $9=="nonsynonymous SNV"' prefix.CLv4_multianno.txt > nonsynonymous_snvs.txt
```


---

## Citation


**Key Software:**
- **BWA:** Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.
- **GATK:** McKenna et al. (2010) The Genome Analysis Toolkit. Genome Research, 20:1297-1303.
- **SAMtools:** Li et al. (2009) The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25:2078-2079.
- **ANNOVAR:** Wang et al. (2010) ANNOVAR: functional annotation of genetic variants. Nucleic Acids Research, 38:e164.

---

**Last Updated:** January 2025
**Version:** 1.0
