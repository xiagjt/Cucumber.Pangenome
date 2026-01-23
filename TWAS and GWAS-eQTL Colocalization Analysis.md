# TWAS and GWAS-eQTL Colocalization Analysis

## Overview
This pipeline performs Transcriptome-Wide Association Study (TWAS) to identify gene expression-trait associations and conducts colocalization analysis to determine if GWAS and eQTL signals share the same causal variant.

## Table of Contents
- [Part 1: TWAS Analysis](#part-1-twas-analysis)
- [Part 2: GWAS-eQTL Colocalization Analysis](#part-2-gwas-eqtl-colocalization-analysis)
- [Software Requirements](#software-requirements)


---

## Part 1: TWAS Analysis

### Overview
Transcriptome-Wide Association Study (TWAS) tests for associations between gene expression levels and phenotypic traits while accounting for population structure using a kinship matrix.

### Workflow

#### Step 1: Prepare Phenotype Data

Extract phenotype values for specific trait from raw data file.

```bash
awk -v colname=FL2019F 'BEGIN {FS="\t"; header=1} \
  NR==header {for (i=1; i<=NF; i++) if ($i==colname) colnum=i} \
  NR>=header {if (colnum) print $1"\t"$colnum}' \
  ../../../../data/Traits.Raw.xls > ./Exp.FL2019F.Trait.txt
```

**Parameters:**
- `colname`: Name of the trait column to extract
- Input: `Traits.Raw.xls` (tab-delimited trait matrix)
- Output: `Exp.FL2019F.Trait.txt` (sample ID + trait value)

**Output format:**
```
SampleID    Trait_Value
Sample001   12.5
Sample002   15.3
Sample003   -9
```

#### Step 2: Filter and Sort Phenotype Data

Select top/bottom samples and remove missing values (-9).

```bash
sort -k 2,2g Exp.FL2019F.Trait.txt | \
  awk '{if($2!="-9") print}' | \
  head -n 50 > Exp.FL2019F.Trait.pole.txt
```

**Purpose:**
- Sort by trait value (ascending)
- Remove missing data (coded as -9)
- Select top 50 samples (can be adjusted)

**Note:** This step performs extreme phenotype sampling, which can increase statistical power for quantitative traits.

---

#### Step 3: Calculate Kinship Matrix

Compute genetic relationship matrix to account for population structure.

```bash
emmax-kin -v -d 10 FL2019F.pole.plink
```

**Parameters:**
- `-v`: Verbose mode (print progress)
- `-d 10`: Use IBS (Identity By State) kinship calculation

**Input:** 
- `FL2019F.pole.plink` (PLINK binary files: .bed, .bim, .fam)

**Output:**
- `FL2019F.pole.plink.BN.kinf` - Balanced kinship matrix
- `FL2019F.pole.plink.hBN.kinf` - Inverse of kinship matrix
- `FL2019F.pole.plink.aIBS.kinf` - IBS kinship matrix

**Kinship Matrix Format:**
```
0.523  0.102  0.089  ...
0.102  0.612  0.134  ...
0.089  0.134  0.598  ...
```

---

#### Step 4: Prepare Gene Expression Data

Order expression matrix to match sample order in phenotype file.

```R
# Read expression matrix (TPM normalized and quantile normalized)
expr_matrix1 <- read.table("../../../data/TPM.qqnorm.xls", 
                          row.names = 1, 
                          sep = "\t", 
                          header = TRUE)

# Read sample IDs from PLINK file
sample_ids <- readLines(file = "FL2019F.pole.plink.id")

# Reorder expression matrix to match phenotype sample order
ordered_matrix1 <- expr_matrix1[sample_ids, , drop = FALSE]

# Write ordered matrix
write.table(ordered_matrix1,
           file = "TPM.qqnorm.ordered.transposed.xls", 
           quote = FALSE, 
           sep = "\t")

# Add header with ID column name
sed -i '1s/^/ID\t/' TPM.qqnorm.ordered.transposed.xls
```

**Input:**
- `TPM.qqnorm.xls` - Gene expression matrix (samples × genes)
  - TPM normalized
  - Quantile normalized (qqnorm)

**Output:**
- `TPM.qqnorm.ordered.transposed.xls` - Ordered expression matrix

**Expression Matrix Format:**
```
ID          Gene001    Gene002    Gene003    ...
Sample001   5.234      12.456     0.123      ...
Sample002   6.789      11.234     0.456      ...
Sample003   4.567      13.678     0.234      ...
```

---

#### Step 5: Run TWAS Analysis

Perform association testing between gene expression and trait using EMMAX.

```R
library(tidyverse)
library(cpgen)

# Load phenotype data
y <- read.table(file = "Exp.FL2019F.Trait.pole.txt", 
               sep = "\t", 
               header = TRUE) 
y <- y[,-1] %>% as.numeric()  # Remove sample ID column

# Load expression matrix
M <- as.matrix(read.table(file = "TPM.qqnorm.ordered.transposed.xls", 
                         sep = "\t", 
                         header = TRUE))

# Transpose and convert to numeric matrix
Mt <- apply(as.matrix(M[,-1]), 1, as.numeric)
rownames(Mt) <- colnames(M[,-1])

# Filter genes: keep genes with TPM >= 1 in at least 5% of samples
Mt <- Mt[apply(Mt >= 1, 1, sum) >= ncol(Mt) * 0.05, ]

# Load kinship matrix
X <- as.matrix(read.table(file = "FL2019F.pole.plink.BN.kinf", 
                         sep = "\t", 
                         header = FALSE))

# Run TWAS with log2 transformation
res <- cGWAS.emmax(log2(y), 
                   log2(t(Mt + 1)),  # Add 1 to avoid log(0)
                   X, 
                   dom = FALSE,      # No dominance effect
                   verbose = TRUE)

# Prepare results
final <- data.frame(gene = rownames(Mt), 
                   beta = res$beta,      # Effect size
                   p_value = res$p_value, # P-value
                   se = res$se)          # Standard error

# Save results
write.table(final,
           file = "./FL2019F.TWAS.xls",
           quote = FALSE,
           row.names = FALSE,
           sep = "\t")
```

**Key Parameters:**
- **Expression filtering**: `TPM >= 1` in `>= 5%` of samples
- **Transformation**: `log2(TPM + 1)` for both trait and expression
- **Model**: EMMAX (Efficient Mixed-Model Association eXpedited)
- **Kinship**: Controls for population structure

**Alternative Transformations:**
```R
# Option 1: No transformation
res <- cGWAS.emmax(y, t(Mt), X, dom=FALSE, verbose=TRUE)

# Option 2: Log10 transformation
res <- cGWAS.emmax(log10(y), log10(t(Mt)), X, dom=FALSE, verbose=TRUE)

# Option 3: Log2 transformation (recommended)
res <- cGWAS.emmax(log2(y), log2(t(Mt+1)), X, dom=FALSE, verbose=TRUE)
```

---

### TWAS Output

**File:** `FL2019F.TWAS.xls`

```
gene              beta           p_value        se
CsaV4_1G000010   0.234567      1.234e-05      0.045678
CsaV4_1G000020   -0.123456     2.345e-04      0.034567
CsaV4_1G000030   0.456789      8.901e-08      0.067890
```

**Columns:**
- **gene**: Gene identifier
- **beta**: Effect size (change in trait per unit change in expression)
- **p_value**: Association p-value
- **se**: Standard error of beta

**Significance Threshold:**
- Bonferroni correction: `0.05 / number_of_genes_tested`
- Example: For 20,000 genes, threshold = `2.5e-06`

---

## Part 2: GWAS-eQTL Colocalization Analysis

### Overview
Colocalization analysis determines whether GWAS and eQTL signals in the same genomic region are driven by the same causal variant, using Bayesian inference (COLOC method).

### Workflow

#### Step 1: Identify Significant GWAS Loci and Perform LD Analysis

Extract genome-wide significant SNPs and identify LD blocks.

```bash
# Run LD analysis on significant GWAS SNPs (p < 3.026268e-05)
/home/Software.Repository/plink/plink \
  --bfile /path/to/genotypes.impute.GWAS \
  --extract SV.GWAS.FL2019S-1.Filter.GWAS.IDs \
  --r2 \
  --ld-window-r2 0.2 \
  --ld-window-kb 100 \
  --threads 100 \
  --memory 20000 \
  --allow-extra-chr \
  --out SV.GWAS.FL2019S-1.Filter.GWAS.IDs.LD
```

**Parameters:**
- `--bfile`: PLINK binary genotype files
- `--extract`: List of significant GWAS SNP IDs
- `--r2`: Output r² values
- `--ld-window-r2 0.2`: Minimum r² threshold (0.2)
- `--ld-window-kb 100`: LD window size (100 kb)

**Input:**
- `SV.GWAS.FL2019S-1.Filter.GWAS.IDs` - List of significant SNP IDs (one per line)

**Output:**
- `.ld` file with LD information

---

#### Step 2: Merge Adjacent LD Regions

Merge loci within 1 Mb to define candidate regions.

```bash
# Convert LD file to BED format (if needed)
# Then merge adjacent regions

/home/Software.Repository/Bedtools/bedtools.static.binary merge \
  -i SV.GWAS.FL2019S-1.Filter.GWAS.IDs.LD.ld.bed \
  -d 1000000 \
  > SV.GWAS.FL2019S-1.Filter.GWAS.IDs.LD.ld.merge.bed

cat SV.GWAS.FL2019S-1.Filter.GWAS.IDs.LD.ld.merge.bed
```

**Parameters:**
- `-i`: Input BED file
- `-d 1000000`: Maximum distance (1 Mb) for merging

**Output:**
```
chr     start       end
3       35818797    37139844
```

**Purpose:** Define candidate regions containing GWAS signals for downstream colocalization analysis.

---

#### Step 3: Extract GWAS Summary Statistics for Candidate Region

Extract all SNPs from GWAS results that fall within the merged LD region.

```bash
# Extract header
head -n 1 ../data/SV.GWAS.FL2019S-1.assoc.txt > INS22773.R.1MB.GWAS.txt

# Extract SNPs in the region (chr3:35818797-37139844)
awk '{if($1 == 3 && $3 >= 35818797 && $3 <= 37139844) print}' \
  ../data/SV.GWAS.FL2019S-1.assoc.txt >> INS22773.R.1MB.GWAS.txt
```

**Output format:**
```
chr  rs         ps        n_miss  allele1  allele0  af     beta          se           logl_H1      l_remle      p_wald
3    DEL28770   35818797  0       G        A        0.459  -1.970e-01    4.281e-02    -9.798e+01   1.478e+00    7.799e-06
3    INS22633   35825251  0       G        A        0.203   9.458e-02    6.539e-02    -1.060e+02   4.188e+00    1.497e-01
```

**Columns:**
- **chr**: Chromosome
- **rs**: SNP/SV identifier
- **ps**: Physical position (bp)
- **n_miss**: Number of missing genotypes
- **allele1/allele0**: Reference/alternative alleles
- **af**: Allele frequency
- **beta**: Effect size
- **se**: Standard error
- **p_wald**: Wald test p-value

---

#### Step 4: Extract Cis-eQTLs for Genes in Candidate Region

Identify genes in the region and extract their cis-eQTL results.

##### 4a. Identify Genes in Region

```bash
# Extract gene IDs from GTF file for chr3:35818797-37139844
awk '{if($1 == "chr3" && $4 >= 35818797 && $5 <= 37139844) print}' \
  ../../../../../LncRNA.dir/data/genome/CLv4.Gene.gtf | \
  awk '{print $12}' | \
  sed -e 's/"//g' -e 's/;//g' | \
  sort | uniq > INS22773.R.1MB.gene.id

cat INS22773.R.1MB.gene.id
```

**Output:**
```
CsaV4_3G003710
CsaV4_3G003720
CsaV4_3G003730
```

##### 4b. Extract eQTL Results

```bash
# For each gene, extract eQTL results for GWAS SNPs
for i in $(cat INS22773.R.1MB.gene.id)
do
  cut -f 2 INS22773.R.1MB.GWAS.txt | \
    grep -Fwf - ../data/Emmax.dir/qqnorm.${i}.txt.ps | \
    awk -v a=${i} '{print $0"\t"a}' - | \
    awk 'NR==FNR{a[$2]=$5; next} $1 in a{print $0"\t"a[$1]}' \
      ../data/genotypes.impute.GWAS.freq_stat.frq - | \
    awk 'NR==FNR{a[$3]=$1"\t"$2; next} $1 in a{print $0"\t"a[$1]}' \
      ../data/SV.xls - >> INS22773.R.1MB.eQTL.txt
done

# Add header
sed -i '1i rs\tbeta\tSE\tp_value\tgene_id\tmaf\tchr\tpos' INS22773.R.1MB.eQTL.txt
```

**Output format:**
```
rs         beta           SE            p_value       gene_id           maf    chr    pos
DEL28770   0.1537600294   0.07848076    0.05161607    CsaV4_3G003710   0.463  chr3   35818797
INS22633   -0.1376663086  0.13837626    0.32112063    CsaV4_3G003710   0.202  chr3   35825251
```

**Columns:**
- **rs**: SNP/SV identifier
- **beta**: eQTL effect size
- **SE**: Standard error
- **p_value**: eQTL p-value
- **gene_id**: Gene identifier
- **maf**: Minor allele frequency
- **chr**: Chromosome
- **pos**: Position

---

#### Step 5: Colocalization Analysis Using COLOC

Test whether GWAS and eQTL signals share the same causal variant.

```R
library(coloc)
library(tidyverse)

# Load eQTL and GWAS data
eqtl <- read.table(file = "INS22773.R.1MB.eQTL.txt", 
                   header = TRUE, 
                   as.is = TRUE) 
gwas <- read.table(file = "INS22773.R.1MB.GWAS.txt", 
                   header = TRUE, 
                   as.is = TRUE)

# Merge datasets by SNP ID
input <- merge(eqtl, gwas, 
              by = "rs", 
              all = FALSE, 
              suffixes = c("_eqtl", "_gwas")) 

# Remove duplicate SNPs (keep best p-value)
input1 <- input[order(input$rs, input$p_value), ]
un_ids <- unique(input1$rs)
index <- match(un_ids, input$rs)
input2 <- input1[index, ]

# Run colocalization analysis
result <- coloc.abf(
  dataset1 = list(
    pvalues = as.numeric(input2$p_wald),  # GWAS p-values
    type = "quant",                       # Quantitative trait
    N = 97,                               # Sample size
    snp = input2$rs
  ), 
  dataset2 = list(
    pvalues = as.numeric(input2$p_value), # eQTL p-values
    type = "quant",                       # Quantitative trait
    N = 97,                               # Sample size
    snp = input2$rs
  ), 
  MAF = as.numeric(input2$maf)           # Minor allele frequency
)

# Extract and format results
need_result <- result$results %>% 
  dplyr::arrange(desc(SNP.PP.H4))

need_result1 <- input2 %>% 
  select(snp = "rs", "gene_id") %>% 
  merge(need_result, by = "snp", all = FALSE) %>% 
  dplyr::arrange(desc(SNP.PP.H4))

head(need_result1)

# Save results
write.table(need_result1,
           file = "INS22773.PP.H4.txt", 
           row.names = FALSE, 
           quote = FALSE, 
           sep = "\t")
```

---

### Colocalization Output

**File:** `INS22773.PP.H4.txt`

```
snp         gene_id          SNP.PP.H0    SNP.PP.H1    SNP.PP.H2    SNP.PP.H3    SNP.PP.H4
DEL28770    CsaV4_3G003710   0.001        0.023        0.045        0.067        0.864
INS22633    CsaV4_3G003710   0.012        0.234        0.123        0.234        0.397
```

**Columns:**
- **snp**: SNP/SV identifier
- **gene_id**: Gene identifier
- **SNP.PP.H0**: Posterior probability - Neither trait has genetic association
- **SNP.PP.H1**: Posterior probability - Only GWAS trait has association
- **SNP.PP.H2**: Posterior probability - Only eQTL trait has association
- **SNP.PP.H3**: Posterior probability - Both traits associated, different causal variants
- **SNP.PP.H4**: Posterior probability - Both traits associated, shared causal variant

---

## Software Requirements

### Required Software

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| EMMAX | - | Kinship matrix calculation | [EMMAX Download](http://csg.sph.umich.edu/kang/emmax/) |
| PLINK | 1.9+ | LD analysis | [PLINK Download](https://www.cog-genomics.org/plink/) |
| BEDTools | 2.0+ | Genomic region operations | [BEDTools Download](https://bedtools.readthedocs.io/) |
| R | 4.0+ | Statistical analysis | [R Project](https://www.r-project.org/) |

### Required R Packages

```R
# Install from CRAN
install.packages(c("tidyverse", "coloc"))

# Install cpgen (for TWAS)
# If not on CRAN, check package repository or GitHub
install.packages("cpgen")
```

---
## Citation

**COLOC method:**
```
Giambartolomei et al. (2014) Bayesian test for colocalisation between pairs 
of genetic association studies using summary statistics. 
PLoS Genetics 10(5): e1004383.
```

**EMMAX:**
```
Kang et al. (2010) Variance component model to account for sample structure 
in genome-wide association studies. Nature Genetics 42:348-354.
```

---
