# Structural Variant Detection and Graph-Pangenome Genotyping Pipeline

## Overview
Comprehensive pipeline for detecting structural variants (SVs) using multiple callers, constructing graph-based pangenomes, and performing SV genotyping across populations. This multi-caller approach maximizes sensitivity and specificity of SV detection.

## Table of Contents
- [Step 0: Read Mapping](#step-0-read-mapping)
- [Step 1: CuteSV](#step-1-cutesv)
- [Step 2: PBSV](#step-2-pbsv)
- [Step 3: Sniffles](#step-3-sniffles)
- [Step 4: SVIM](#step-4-svim)
- [Step 5: MUMmer-Syri (Assembly-based)](#step-5-mummer-syri-assembly-based)
- [Step 6: SV Merging](#step-6-sv-merging)
- [Step 7: Graph-Pangenome Construction](#step-7-graph-pangenome-construction)
- [Step 8: Graph-Pangenome Genotyping](#step-8-graph-pangenome-genotyping)
- [Software Requirements](#software-requirements)


---

## Step 0: Read Mapping

### Overview
Map PacBio long reads to reference genome CLv4.0. Two mapping approaches are used:
- **pbmm2** for PBSV caller (optimized for PacBio data)
- **minimap2** for CuteSV, Sniffles, and SVIM callers

### 0.1 Mapping for PBSV (using pbmm2)

```bash
# Map HiFi reads using pbmm2
pbmm2 align \
  --preset HiFi \
  --sort \
  -j 80 \
  --log-level INFO \
  CLv4.0.fa \
  CG37.hifi_reads.bam \
  CG37.pbmm2.bam
```

**Parameters:**
- `--preset HiFi`: Optimized for PacBio HiFi reads
- `--sort`: Output sorted BAM
- `-j 80`: Number of threads
- `--log-level INFO`: Logging verbosity

**For PacBio CLR reads (from previous study):**
```bash
pbmm2 align \
  --preset CCS \
  --sort \
  -j 80 \
  CLv4.0.fa \
  sample.clr_reads.bam \
  sample.pbmm2.bam
```

### 0.2 Mapping for Other Callers (using minimap2)

```bash
# Map HiFi reads to reference
minimap2 -ax map-hifi -t 80 \
  CLv4.0.fa \
  CG37.hifi_reads.bam \
  > CG37.minimap2.raw.sam

# For PacBio CLR reads (from previous study)
minimap2 -ax map-pb -t 80 \
  CLv4.0.fa \
  sample.clr_reads.bam \
  > sample.minimap2.raw.sam

# Convert SAM to BAM
samtools view -bh -@ 20 \
  CG37.minimap2.raw.sam \
  > CG37.minimap2.raw.bam

# Sort BAM file
samtools sort -@ 20 \
  -T CG37.minimap2.raw.bam \
  -o CG37.minimap2.sort.bam \
  CG37.minimap2.raw.bam

# Index sorted BAM
samtools index -@ 20 CG37.minimap2.sort.bam
```

### Parameters Explained

| Parameter | Tool | Description |
|-----------|------|-------------|
| `-ax map-hifi` | minimap2 | Preset for PacBio HiFi reads (Q20+, high accuracy) |
| `-ax map-pb` | minimap2 | Preset for PacBio CLR reads (older technology) |
| `-t 80` | minimap2 | Number of threads |
| `-j 80` | pbmm2 | Number of threads |
| `--preset HiFi` | pbmm2 | Optimized for HiFi reads |
| `--preset CCS` | pbmm2 | For circular consensus sequencing (CLR) |
| `-bh` | samtools view | Output BAM format with header |
| `-@ 20` | samtools | Number of threads |

### Input Files
- `CLv4.0.fa` - Reference genome (near-complete cucumber genome assembly)
- `*.hifi_reads.bam` - PacBio HiFi reads (125 accessions)
- `*.clr_reads.bam` - PacBio CLR reads (8 accessions from previous study)

### Output Files

**For PBSV:**
- `*.pbmm2.bam` - Sorted, aligned BAM file
- `*.pbmm2.bam.bai` - BAM index file

**For other callers:**
- `*.minimap2.sort.bam` - Sorted BAM file
- `*.minimap2.sort.bam.bai` - BAM index file

### Quality Checks
```bash
# Check mapping statistics
samtools flagstat sample.minimap2.sort.bam

# Check coverage
samtools depth sample.minimap2.sort.bam | \
  awk '{sum+=$3} END {print "Average coverage:", sum/NR}'
```

**Expected Results:**
- **HiFi reads:** Mapping rate >95%, average coverage 20-40×
- **CLR reads:** Mapping rate >90%, average coverage 30-60×

### Data Summary
- **Total accessions:** 133 (125 HiFi + 8 CLR)
- **HiFi data:** 125 newly sequenced accessions
- **CLR data:** 8 accessions from previous study (published data)

---

## Step 1: CuteSV

### Overview
Long-read-based SV caller optimized for PacBio HiFi and ONT reads. Good for detecting insertions and deletions.

### 1.1 Calling

```bash
cuteSV \
  CG104.minimap2.sort.bam \
  NearC.V1.fa \
  CG104.cutesv.raw.vcf \
  CG104.out \
  --max_cluster_bias_INS 1000 \
  --diff_ratio_merging_INS 0.9 \
  --max_cluster_bias_DEL 1000 \
  --diff_ratio_merging_DEL 0.5 \
  --threads 80
```

**Parameters:**
- `--max_cluster_bias_INS 1000`: Maximum distance for clustering insertion signals (bp)
- `--diff_ratio_merging_INS 0.9`: Minimum similarity for merging insertions (0-1)
- `--max_cluster_bias_DEL 1000`: Maximum distance for clustering deletion signals (bp)
- `--diff_ratio_merging_DEL 0.5`: Minimum similarity for merging deletions (0-1)
- `--threads 80`: Number of CPU threads

### 1.2 Filtering

```bash
for f in *cutesv.raw.vcf
do
  cat $f | perl -ane '
    if(/^#/){
      print $_;
    } else {
      my $len;
      if(/SVLEN=(\S+?)\;/){$len=$1;}
      my $re;
      if(/RE=(\S+?)\;/){$re=$1;}
      if($re>=4 && abs($len)<20000 && $_!~/BND/ && /PASS/){
        print $_;
      }
    }' - > ${f}.filter.vcf
done
```

**Filtering Criteria:**
- **RE ≥ 4**: Read support ≥ 4 (number of reads supporting SV)
- **|SVLEN| < 20,000**: SV length < 20 kb (shorter than 20,000 bp)
- **No BND**: Exclude breakends (complex SVs)
- **PASS**: Only variants passing all filters

**Output:**
- `*.cutesv.raw.vcf.filter.vcf` - Filtered SV calls

---

## Step 2: PBSV

### Overview
Official PacBio structural variant caller, optimized for PacBio data. Excellent for complex SVs and inversions.

### 2.1 Calling

```bash
# First, discover SV signatures (if not already done)
pbmm2 align --preset HiFi --sort NearC.V1.fa R70_2.hifi.bam R70_2.pbmm2.bam
pbsv discover R70_2.pbmm2.bam R70_2.pbmm2.bam.svsig.gz

pbmm2 align --preset HiFi --sort NearC.V1.fa R70.hifi.bam R70.pbmm2.bam
pbsv discover R70.pbmm2.bam R70.pbmm2.bam.svsig.gz

# Call SVs
pbsv call -j 100 --hifi \
  ../Bams.dir/NearC.V1.fa \
  R70_2.pbmm2.bam.svsig.gz R70.pbmm2.bam.svsig.gz \
  R70.vcf
```

**Parameters:**
- `-j 100`: Number of threads
- `--hifi`: Optimize for HiFi reads

### 2.2 Filtering

```bash
for f in *.vcf
do
  bcftools view -i 'DP>=4' $f | perl -ane '
    if(/^#/){
      print $_;
    } else {
      s/INS.DUP/DUP/mg;
      if($_!~/PASS/){next;}
      if(/BND/){next;}
      my $len;
      if(/SVLEN=(\S+?)\;/){$len=$1;}
      if(abs($len)<20000){
        print $_;
      }
    }' - > ${f}.filter.vcf
done
```

**Filtering Criteria:**
- **DP ≥ 4**: Depth ≥ 4
- **PASS**: Only passing variants
- **No BND**: Exclude breakends
- **|SVLEN| < 20,000**: SV length < 20 kb
- **Rename**: Convert INS.DUP to DUP for consistency

---

## Step 3: Sniffles

### Overview
Fast and sensitive SV caller for long reads. Good balance between speed and accuracy.

### 3.1 Calling

```bash
sniffles \
  --input CG104.minimap2.sort.bam \
  --vcf CG104.sniffles.raw.vcf \
  --reference NearC.V1.fa \
  -t 10
```

**Parameters:**
- `--input`: Input BAM file
- `--vcf`: Output VCF file
- `--reference`: Reference genome (optional but recommended)
- `-t`: Number of threads

### 3.2 Filtering

```bash
for f in *sniffles.raw.vcf
do
  bcftools view -i 'QUAL>=10' $f | perl -ane '
    if(/^#/){
      print $_;
    } else {
      my $len;
      if(/SVLEN=(\S+?)\;/){$len=$1;}
      my @g = split /\:/, $_;
      my $DP = $g[-1]+$g[-2];
      if($DP>=4 && abs($len)<20000 && $_!~/BND/ && /PASS/){
        print $_;
      }
    }' - > ${f}.filter.vcf
done
```

**Filtering Criteria:**
- **QUAL ≥ 10**: Quality score ≥ 10
- **DP ≥ 4**: Read depth ≥ 4 (calculated from genotype fields)
- **|SVLEN| < 20,000**: SV length < 20 kb
- **No BND**: Exclude breakends
- **PASS**: Only passing variants

---

## Step 4: SVIM

### Overview
SV caller based on alignment signatures. Good for all types of SVs including inversions.

### 4.1 Calling

```bash
svim alignment \
  --sample CG104 \
  CG104.out \
  CG104.minimap2.sort.bam \
  NearC.V1.fa
```

**Parameters:**
- `--sample`: Sample name
- Output directory: `CG104.out`
- Input: BAM file and reference genome

### 4.2 Filtering

```bash
for f in *out/variants.vcf
do
  bcftools view -i 'QUAL>=10 && DP>=4' $f | perl -ane '
    if(/^#/){
      print $_;
    } else {
      if($_=~/BND/ || $_!~/PASS/){next}
      my $len;
      if($_!~/INV/i){
        if(/SVLEN=(\S+?)\;/){$len=$1;}
        if(abs($len)<20000){print $_;}
      } else {
        if($_=~/BND/ || $_!~/PASS/){next}
        my $end;
        my @g = split /\t/, $_;
        if(/END=(\d+)\;/){$end=$1;}
        if(abs($end-$g[1])<20000){print $_;}
      }
    }' - > ${f}.filter.vcf
done
```

**Filtering Criteria:**
- **QUAL ≥ 10**: Quality score ≥ 10
- **DP ≥ 4**: Depth ≥ 4
- **PASS**: Only passing variants
- **No BND**: Exclude breakends
- **For non-inversions**: |SVLEN| < 20,000
- **For inversions**: |END - POS| < 20,000

---

## Step 5: MUMmer-Syri (Assembly-based)

### Overview
Genome alignment-based SV detection. Gold standard for comparing genome assemblies.

### 5.1 Genome Alignment with MUMmer

```bash
nucmer \
  --mum \
  -c 100 \
  -l 50 \
  -g 500 \
  -t 40 \
  -p CG25.mummer \
  CLv4.0.chr.fa \
  CG25.genome.chr.fa
```

**Parameters:**
- `--mum`: Use maximal unique matches (only unique matches in both genomes)
- `-c 100`: Minimum cluster length (100 bp)
- `-l 50`: Minimum match length (50 bp)
- `-g 500`: Maximum gap between matches (500 bp)
- `-t 40`: Number of threads
- `-p`: Output prefix

**Purpose:**
These parameters are optimized for detecting structural variations between closely related genomes:
- `-c 100` and `-l 50` ensure sufficient alignment quality
- `-g 500` allows detection of larger indels and rearrangements
- `--mum` provides high specificity by requiring unique matches

**Output:**
- `CG25.mummer.delta` - Alignment delta file containing all alignments

### 5.2 SV Detection with Syri

```bash
syri \
  --nc 40 \
  -c CG28.mummer.coords \
  -d CG28.mummer.filter.delta \
  -r NearC.V1.chr.fa \
  -q CG28.genome.chr.fa \
  --prefix CG28.syri.mummer
```

**Parameters:**
- `--nc 40`: Number of CPU cores
- `-c`: Coordinates file from delta-filter/show-coords
- `-d`: Delta file from nucmer
- `-r`: Reference genome
- `-q`: Query genome
- `--prefix`: Output prefix

**Note:** May need to run `delta-filter` and `show-coords` on delta file first:
```bash
delta-filter -1 CG28.mummer.delta > CG28.mummer.filter.delta
show-coords -T -H CG28.mummer.filter.delta > CG28.mummer.coords
```

---

## Step 6: SV Merging

### Overview
Merge SV calls from multiple callers (5 callers total: CuteSV, PBSV, Sniffles, SVIM, and Syri) using SURVIVOR to create a consensus callset. **Important**: Only SVs < 20 kb are merged; SVs ≥ 20 kb detected by MUMmer are preserved separately.

### 6.1 Merge SVs < 20 kb from Five Callers

```bash
/bioinfo/guanjt/Software.dir/SURVIVOR-master/Debug/SURVIVOR merge \
  file.CG16.txt \
  1000 \
  2 \
  1 \
  1 \
  0 \
  20 \
  CG16.merge.vcf
```

### Parameters

| Parameter | Position | Value | Description |
|-----------|----------|-------|-------------|
| Input file list | 1 | `file.CG16.txt` | Text file with VCF file paths (one per line) |
| Max distance | 2 | `1000` | Maximum distance between breakpoints (1 kb) |
| Min callers | 3 | `2` | Minimum 2 callers required to support SV |
| Type match | 4 | `1` | Same SV type required (DEL, INS, etc.) |
| Strand match | 5 | `1` | Same strand required |
| Distance estimate | 6 | `0` | Use average distance |
| Min size | 7 | `20` | Minimum SV size (20 bp) |
| Output | 8 | `CG16.merge.vcf` | Output merged VCF file |

### Input File Format

**file.CG16.txt** (list of 5 filtered VCF files):
```
CG16.cutesv.raw.vcf.filter.vcf
CG16.pbsv.vcf.filter.vcf
CG16.sniffles.raw.vcf.filter.vcf
CG16.out/variants.vcf.filter.vcf
CG16.syri.mummer.vcf
```

**Note:** 
- First 4 files: From long-read based callers (CuteSV, PBSV, Sniffles, SVIM)
- Last file: From MUMmer-Syri (assembly-based caller)
- All input VCFs should contain only SVs < 20 kb (filtered in previous steps)

### 6.2 Preserve Large SVs from MUMmer

```bash
# Extract SVs >= 20 kb from MUMmer-based Syri results
bcftools view CG16.syri.mummer.vcf | \
  awk 'BEGIN {FS="\t"; OFS="\t"} 
    /^#/ {print; next}
    {
      if(/SVLEN=([^;]+)/) {
        len = substr($0, RSTART+6, RLENGTH-6);
        if(len ~ /^-/) len = -len;
        if(len >= 20000) print;
      }
    }' > CG16.large_svs.vcf
```

### 6.3 Combine Final SV Callset

```bash
# Combine merged SVs (<20kb) with large SVs (>=20kb)
bcftools concat \
  CG16.merge.vcf \
  CG16.large_svs.vcf \
  -a \
  -o CG16.final.vcf

# Sort and remove duplicates
bcftools sort CG16.final.vcf | \
  bcftools norm -d all -o CG16.final.sorted.vcf
```

### Merging Strategy Details

**For SVs < 20 kb (133 accessions):**
- **Data source:** 5 callers (CuteSV, PBSV, Sniffles, SVIM, Syri from reads/BAM)
- **Max distance:** 1000 bp (1 kb)
- **Min support:** ≥2 callers
- **Type matching:** Required (must agree on DEL, INS, etc.)
- **Rationale:** Reduces false positives while maintaining sensitivity

**For SVs ≥ 20 kb:**
- **Data source:** MUMmer whole-genome alignment
- **No merging:** Preserved as-is from genome assemblies
- **Rationale:** Assembly-based calling is more accurate for large SVs

### Why This Two-Tier Approach?

1. **Read-based callers** (CuteSV, PBSV, Sniffles, SVIM):
   - Excellent for small to medium SVs (50 bp - 20 kb)
   - Multiple callers reduce false positives
   - Limited accuracy for very large SVs

2. **Assembly-based caller** (MUMmer-Syri):
   - Gold standard for large SVs (≥ 20 kb)
   - No merging needed (high confidence)
   - Can detect complex rearrangements

### Quality Metrics

**Check caller support distribution:**
```bash
# Count SVs by number of supporting callers
bcftools query -f '%SUPP\n' CG16.merge.vcf | sort | uniq -c

# Expected output:
#   1234  2    (supported by 2 callers)
#   5678  3    (supported by 3 callers)
#   2345  4    (supported by 4 callers)
#    456  5    (supported by all 5 callers)
```

**Check SV type distribution:**
```bash
bcftools query -f '%SVTYPE\n' CG16.final.sorted.vcf | sort | uniq -c
```

### Alternative Merging Strategies

**Conservative (higher specificity):**
```bash
SURVIVOR merge file.txt 500 3 1 1 0 20 output.vcf
# Require 3+ callers, smaller distance window
```

**Sensitive (higher sensitivity):**
```bash
SURVIVOR merge file.txt 2000 2 1 1 0 20 output.vcf  
# Allow 2 kb distance, still require 2+ callers
```

**Our approach (balanced):**
```bash
SURVIVOR merge file.txt 1000 2 1 1 0 20 output.vcf
# 1 kb distance, 2+ callers - good balance
```

### Output Files

| File | Description | SV Count |
|------|-------------|----------|
| `*.merge.vcf` | Merged SVs < 20 kb from 5 callers | 10,000-15,000 |
| `*.large_svs.vcf` | SVs ≥ 20 kb from MUMmer | 500-1,000 |
| `*.final.sorted.vcf` | Final complete callset | 10,500-16,000 |

### Summary Statistics

**Expected per genome:**
- Total SVs: 10,000-16,000
- SVs < 20 kb (merged): ~95%
- SVs ≥ 20 kb (preserved): ~5%
- Caller support breakdown:
  - 2 callers: ~40-50%
  - 3 callers: ~30-40%
  - 4 callers: ~10-20%
  - 5 callers: ~5-10%

---

## Step 7: Graph-Pangenome Construction

### Overview
Construct variation graph from reference genome and SVs using vg toolkit. Focuses on presence-absence variations (PAVs) and inversions.

### Commands

```bash
threads=80
vcf="variants.vcf.gz"  # Must be bgzipped and indexed
ref="reference.fa"

# 1. Construct graph from reference + VCF
vg construct \
  -t $threads \
  -a \
  -f \
  -S \
  -v $vcf \
  -r $ref \
  > graph.vg

# 2. Build XG index
vg index \
  -t $threads \
  -L \
  -b tmp \
  -x graph.xg \
  graph.vg

# 3. Build GBWT index (haplotype paths)
vg gbwt \
  -d tmp \
  -g graph.gg \
  -x graph.xg \
  -o graph.gbwt \
  -P

# 4. Extract snarls (variation sites)
vg snarls \
  -t $threads \
  --include-trivial \
  graph.xg \
  > graph.snarls

# 5. Build distance index
vg index \
  -b tmp \
  -t $threads \
  -j graph.dist \
  -s graph.snarls \
  graph.vg

# 6. Build minimizer index
vg minimizer \
  -t $threads \
  -i graph.min \
  -g graph.gbwt \
  -d graph.dist \
  graph.xg

# 7. Re-extract snarls (for calling)
vg snarls \
  -t $threads \
  graph.xg \
  > graph.snarls
```

### Parameters Explained

#### vg construct
- `-t`: Number of threads
- `-a`: Index alt paths (important for variants)
- `-f`: Include all variants regardless of phase
- `-S`: Store paths for all samples
- `-v`: Input VCF file (must be bgzipped and indexed)
- `-r`: Reference FASTA file

#### vg index
- `-L`: Build hash-based index
- `-b tmp`: Temporary directory
- `-x`: XG index output
- `-j`: Distance index output
- `-s`: Snarls file (variation sites)

#### vg gbwt
- `-d tmp`: Temporary directory
- `-g`: GG graph output
- `-x`: XG index input
- `-o`: GBWT output (haplotype paths)
- `-P`: Extract sample paths from VCF

#### vg snarls
- `--include-trivial`: Include simple bubbles

#### vg minimizer
- `-i`: Minimizer index output
- `-g`: GBWT input
- `-d`: Distance index input

### Input Requirements

**VCF File Preparation:**
```bash
# Must be bgzipped and indexed
bgzip -c variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz
```

**VCF Requirements:**
- Only PAVs (presence-absence variations) and inversions
- Other SV types (duplications) handled separately
- Properly formatted with GT fields

### Output Files

| File | Description | Size |
|------|-------------|------|
| `graph.vg` | Variation graph | Large |
| `graph.xg` | XG index (succinct graph) | Medium |
| `graph.gbwt` | GBWT index (haplotypes) | Medium |
| `graph.gg` | GG graph (for GBWT) | Small |
| `graph.snarls` | Snarls (variation sites) | Small |
| `graph.dist` | Distance index | Medium |
| `graph.min` | Minimizer index | Large |

---

## Step 8: Graph-Pangenome Genotyping

### Overview
Genotype SVs across populations using graph-based approach. Two methods: vg for PAVs/inversions, Paragraph for duplications.

### 8.1 PAVs and Inversions (using vg)

```bash
threads=80
sample="Sample001"

# 1. Map reads to graph using giraffe
vg giraffe \
  -t $threads \
  -f ${sample}_1.clean.fastq.gz \
  -f ${sample}_2.clean.fastq.gz \
  -x graph.xg \
  -g graph.gg \
  -m graph.min \
  -d graph.dist \
  -H graph.gbwt \
  > ${sample}.gam

# 2. Compute read support (pileup)
vg pack \
  -t $threads \
  -Q 5 \
  -x graph.xg \
  -g ${sample}.gam \
  -o ${sample}.pack

# 3. Call variants
vg call \
  -t $threads \
  -r graph.snarls \
  -k ${sample}.pack \
  -s ${sample} \
  graph.xg \
  > ${sample}.vcf
```

### Parameters Explained

#### vg giraffe (fast read mapping)
- `-t`: Number of threads
- `-f`: FASTQ input files (can specify twice for paired-end)
- `-x`: XG graph index
- `-g`: GG graph
- `-m`: Minimizer index
- `-d`: Distance index
- `-H`: GBWT haplotype index

**Output:** `.gam` file (graph alignment format)

#### vg pack (compute coverage)
- `-t`: Number of threads
- `-Q 5`: Minimum mapping quality (5)
- `-x`: XG index
- `-g`: Input GAM file
- `-o`: Output pack file

**Output:** `.pack` file (read coverage at each node)

#### vg call (variant calling)
- `-t`: Number of threads
- `-r`: Snarls file (variation sites)
- `-k`: Pack file (read support)
- `-s`: Sample name
- Input: XG index

**Output:** `.vcf` file with genotypes

---

### 8.2 Duplications (using Paragraph)

```bash
# Genotype duplications using Paragraph
multigrmpy.py \
  -i $vcf \
  -M 2000 \
  -m $bam \
  -r $ref \
  -o Sample.dir \
  --threads 100
```

### Parameters

- `-i`: Input VCF with duplication SVs
- `-M 2000`: Maximum SV size to genotype (bp)
- `-m`: Input BAM file (aligned reads)
- `-r`: Reference genome FASTA
- `-o`: Output directory
- `--threads`: Number of threads

### Why Separate Methods?

**vg (Graph-based):**
- ✅ Excellent for PAVs and inversions
- ✅ Handles complex variation
- ✅ Fast with giraffe mapper
- ❌ Not optimal for duplications

**Paragraph (Read-depth based):**
- ✅ Designed for copy number variation
- ✅ Accurate for duplications
- ✅ Uses read-depth and pair information
- ❌ Slower than graph-based methods

---

## Software Requirements

### Required Software

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| minimap2 | 2.17+ | Read mapping | [GitHub](https://github.com/lh3/minimap2) |
| samtools | 1.12+ | BAM processing | [samtools.org](http://www.htslib.org/) |
| CuteSV | 2.0+ | SV calling | [GitHub](https://github.com/tjiangHIT/cuteSV) |
| pbsv | 2.9.0+ | PacBio SV calling | [GitHub](https://github.com/PacificBiosciences/pbsv) |
| Sniffles | 2.0+ | SV calling | [GitHub](https://github.com/fritzsedlazeck/Sniffles) |
| SVIM | 2.0+ | SV calling | [GitHub](https://github.com/eldariont/svim) |
| Syri | 1.6+ | SV calling from alignments | [GitHub](https://github.com/schneebergerlab/syri) |
| MUMmer | 4.0+ | Genome alignment | [GitHub](https://github.com/mummer4/mummer) |
| SURVIVOR | 1.0.7+ | SV merging | [GitHub](https://github.com/fritzsedlazeck/SURVIVOR) |
| vg toolkit | 1.40+ | Graph-pangenome | [GitHub](https://github.com/vgteam/vg) |
| Paragraph | 2.3+ | Duplication genotyping | [GitHub](https://github.com/Illumina/paragraph) |
| bcftools | 1.12+ | VCF processing | [GitHub](https://github.com/samtools/bcftools) |

### Installation Tips

**Using conda (recommended):**
```bash
# Create environment
conda create -n sv_pipeline

# Activate environment
conda activate sv_pipeline

# Install tools
conda install -c bioconda minimap2 samtools cutesv sniffles svim mummer survivor bcftools

# Install vg (may need manual compilation)
conda install -c bioconda vg

# PBSV requires PacBio tools
conda install -c bioconda pbsv pbmm2
```

**Manual Installation:**
- Follow each tool's GitHub repository instructions
- Ensure all tools are in PATH
- Test with `--version` or `--help`

---

## Best Practices

### 1. Multi-Caller Strategy

**Callers Used in This Study:**

| Caller | Type | Strength | Priority |
|--------|------|----------|----------|
| **CuteSV** | Long-read | Fast, sensitive for INS/DEL | High |
| **PBSV** | Long-read | Official PacBio tool, optimized | High |
| **Sniffles** | Long-read | Balanced speed/accuracy | High |
| **SVIM** | Long-read | Good for inversions | High |
| **MUMmer-Syri** | Assembly-based | Gold standard, especially for large SVs | High |

**Total: 5 Callers**
- 4 long-read based callers
- 1 assembly-based caller

**Why This Combination?**
- **Long-read callers (4)**: Detect SVs from read alignments
  - Each has different algorithms and strengths
  - Consensus approach reduces false positives
  - Good for SVs 50 bp - 20 kb
  
- **Assembly-based caller (1)**: Detects SVs from genome comparisons
  - Most accurate method, especially for large SVs (≥ 20 kb)
  - Can detect complex rearrangements
  - Serves as gold standard

**Merging Strategy:**
- SVs < 20 kb: Require ≥2/5 callers (reduces false positives)
- SVs ≥ 20 kb: Use MUMmer exclusively (highest accuracy)

**Alternative Combinations (Not Used Here):**

| Read Type | Alternative Callers | Notes |
|-----------|-------------------|-------|
| PacBio HiFi | + Picky, DeBreak | Additional sensitivity |
| Oxford Nanopore | + NanoSV, SVIM-asm | ONT-specific tools |
| Short reads | Manta, DELLY, Lumpy | Different technology |

---

### 2. Filtering Strategy

**Standard Filtering (Used in Study):**
```
- Read support (RE/DP) ≥ 4
- Quality score ≥ 10 (for Sniffles and SVIM)
- PASS filter only
- SV length < 20,000 bp
- Present in ≥2 callers (for merged callset)
```

**Conservative Filtering (Higher Specificity):**
```
- Read support (RE/DP) ≥ 10
- Quality score ≥ 20
- PASS filter only
- Present in ≥3 callers
```

**Sensitive Filtering (Higher Sensitivity):**
```
- Read support (RE/DP) ≥ 3
- Quality score ≥ 5
- Include IMPRECISE calls
- Present in ≥1 caller
```

**Note:** The study used **read depth ≥ 4** as a balanced threshold that:
- Reduces false positives from low-coverage regions
- Maintains high sensitivity for true SVs
- Works well across variable sequencing depths (20-40×)

---

### 3. Size Cutoffs

**Recommended Size Ranges:**

| SV Type | Minimum | Maximum | Rationale |
|---------|---------|---------|-----------|
| Small SVs | 50 bp | 1 kb | High confidence |
| Medium SVs | 1 kb | 10 kb | Most biological SVs |
| Large SVs | 10 kb | 20 kb | May need assembly validation |
| Very Large | >20 kb | - | Use assembly-based methods |

**Our Pipeline Uses:**
- Maximum: 20 kb (good balance)
- Minimum: Implicit (usually 50 bp by caller default)

---

### 4. Performance Optimization

**Thread Allocation:**
```bash
# For servers with 100+ cores:
Mapping: 80 threads
SV calling: 40-80 threads
Merging: Single-threaded (fast anyway)
Graph construction: 80 threads
Genotyping: 40-80 threads
```

**Memory Requirements:**
```
Mapping: 8-16 GB
SV calling: 16-32 GB
Graph construction: 32-64 GB (depends on graph size)
Genotyping: 16-32 GB per sample
```

**Disk Space:**
```
Raw reads: 10-20 GB per sample (HiFi)
BAM files: 15-30 GB per sample
VCF files: 100 MB - 2 GB per sample
Graph files: 5-50 GB (depends on variant count)
```

---

### 5. Quality Control

**Check Mapping Quality:**
```bash
samtools flagstat sample.bam
samtools stats sample.bam | grep ^SN
```

**Check SV Counts:**
```bash
# Count SVs by type
bcftools query -f '%SVTYPE\n' sample.vcf | sort | uniq -c

# Count by size bin
bcftools query -f '%SVLEN\n' sample.vcf | awk '{
  len=($1<0)?-$1:$1; 
  if(len<100) bin="<100"; 
  else if(len<1000) bin="100-1k"; 
  else if(len<10000) bin="1k-10k"; 
  else bin=">10k"; 
  count[bin]++
} END {
  for(b in count) print b, count[b]
}'
```

**Check Caller Overlap:**
```bash
# After merging, check SUPP_VEC field
bcftools query -f '%SUPP_VEC\n' merged.vcf | sort | uniq -c
```

---




## Pipeline Summary

### Recommended Workflow

```
1. Map reads → BAM
   - pbmm2 for PBSV
   - minimap2 for other callers
   ↓
2. Call SVs with 4 long-read callers
   - CuteSV
   - PBSV  
   - Sniffles
   - SVIM
   ↓
3. Call SVs from genome assemblies
   - MUMmer-Syri (5th caller)
   ↓
4. Filter each callset
   - Read depth ≥ 4
   - Length < 20 kb
   - PASS flag
   ↓
5. Merge with SURVIVOR
   - SVs < 20 kb: ≥2 callers
   - SVs ≥ 20 kb: Preserve from MUMmer
   ↓
6. Build graph-pangenome (vg)
   - PAVs and inversions only
   ↓
7. Genotype across population
   - PAVs/inversions: vg giraffe + call
   - Duplications: Paragraph
   ↓
8. Filter and validate
```

### Expected Runtime

**For one sample (30× coverage):**
```
Mapping: 2-4 hours (80 threads)
SV calling (all callers): 4-8 hours total
Filtering: <30 minutes
Merging: <5 minutes
```

**For graph construction (100 samples):**
```
VCF preparation: 1-2 hours
Graph construction: 2-6 hours (80 threads)
Index building: 1-3 hours
```

**For genotyping (per sample):**
```
Read mapping (giraffe): 30-60 minutes
Packing: 10-20 minutes
Calling: 10-20 minutes
Total per sample: 1-2 hours
```

---

## Output Summary

### Per-Sample Outputs

| File | Description | Size |
|------|-------------|------|
| `*.minimap2.sort.bam` | Aligned reads | 15-30 GB |
| `*.cutesv.raw.vcf.filter.vcf` | CuteSV calls | 10-50 MB |
| `*.vcf.filter.vcf` | PBSV calls | 10-50 MB |
| `*.sniffles.raw.vcf.filter.vcf` | Sniffles calls | 10-50 MB |
| `*.out/variants.vcf.filter.vcf` | SVIM calls | 10-50 MB |
| `*.merge.vcf` | Merged calls | 10-50 MB |

### Graph Files

| File | Description | Size |
|------|-------------|------|
| `graph.vg` | Variation graph | 5-20 GB |
| `graph.xg` | XG index | 2-10 GB |
| `graph.gbwt` | Haplotype paths | 1-5 GB |
| `graph.min` | Minimizer index | 3-15 GB |

### Genotyping Outputs

| File | Description | Size |
|------|-------------|------|
| `*.gam` | Graph alignments | 10-30 GB |
| `*.pack` | Read support | 100 MB - 1 GB |
| `*.vcf` | Genotypes | 50-200 MB |

---

## Citation



**Key Software Citations:**

- **minimap2:** Li, H. (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100.
- **CuteSV:** Jiang et al. (2020) Long-read-based human genomic structural variation detection with cuteSV. Genome Biology, 21:189.
- **Sniffles:** Sedlazeck et al. (2018) Accurate detection of complex structural variations using single-molecule sequencing. Nature Methods, 15:461-468.
- **SURVIVOR:** Jeffares et al. (2017) Transient structural variations have strong effects on quantitative traits and reproductive isolation in fission yeast. Nature Communications, 8:14061.
- **vg:** Garrison et al. (2018) Variation graph toolkit improves read mapping by representing genetic variation in the reference. Nature Biotechnology, 36:875-879.

---

---

**Last Updated:** January 2025
**Version:** 1.0
