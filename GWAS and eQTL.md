This section describes the genetic mapping workflow, including **Genome-Wide Association Studies (GWAS)** and **expression Quantitative Trait Loci (eQTL)** analysis. These methods are used to identify genetic variants (SNPs/SVs) associated with physical traits and gene expression levels.

---

# Genetic Mapping and eQTL Analysis Pipeline

## 1. Environment Setup

We use **GEMMA** for standard GWAS and **EMMAX** for high-efficiency eQTL mapping.

```bash
# Create environment for genetic mapping
conda create -n mapping_env -c bioconda gemma emmax samtools -y

# Activate the environment
conda activate mapping_env

```

---

## 2. Genome-Wide Association Study (GWAS)

Using **GEMMA** (Genome-wide Efficient Mixed Model Association) to perform association mapping with a Linear Mixed Model (LMM) to account for population structure.

```bash
###### Run GWAS using GEMMA
gemma -bfile Variants.Genome.Final422.Filter \
      -lmm 1 \
      -k ./Related.Matrix/kin.sXX.txt \
      -n 27 \
      -c ./PCA.dir/PCA.Top10.GEMMA.xls \
      -o Raw.GEMMA.GWAS.TC_A

```

### Parameter Description:

* `-bfile`: Prefix for PLINK binary files (`.bed`, `.bim`, `.fam`).
* `-lmm 1`: Specifies the Wald test for the Linear Mixed Model.
* `-k`: Relatedness matrix (Kinship) to control for population stratification.
* `-n 27`: Specifies the column index of the phenotype in the `.fam` file.
* `-c`: Covariate file (e.g., Top 10 PCs from PCA) to control for population structure.
* `-o`: Prefix for the output files.

---

## 3. expression Quantitative Trait Loci (eQTL) Analysis

eQTL analysis identifies genetic variants that regulate the expression levels of genes. We use **EMMAX** for its efficiency in handling thousands of expression phenotypes.

### 3.1 Bulk Processing Script

The following script generates a batch execution file (`Run.sh`) to run EMMAX across multiple gene expression phenotypes.

```bash
# Set paths for kinship and genotype data
kinship="Kinship.dir/Merge.SV.SNP.MAF.aIBS.kinf"
tped="Kinship.dir/Merge.SV.SNP.MAF"

# Loop through all phenotype files (gene expression matrices)
for f in Phenos.dir/*.txt
do
    echo "emmax-intel64 -v -d 10 -t $tped -p $f -c PCA.Top10.xls -k $kinship -o $f" >> Run.sh
done

# Execute the batch jobs
chmod +x Run.sh
./Run.sh

```

### Parameter Description:

* `-v`: Verbose mode.
* `-d 10`: Precision for the EM algorithm.
* `-t`: Prefix for the genotype files in `tped/tfam` format.
* `-p`: Phenotype file (individual gene expression levels).
* `-k`: Kinship matrix.
* `-c`: Covariates (PCA) to reduce false positives.

---
