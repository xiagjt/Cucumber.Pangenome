# Hi-C Data Processing Pipeline

## 📋 Pipeline Overview

This pipeline processes Hi-C data from raw reads to normalized contact maps. It includes read alignment, contact matrix generation, ICE normalization, and visualization. The workflow uses **HiC-Pro** for initial processing and **HiCExplorer** for downstream analysis.

---

## 🔧 Requirements

### Software Dependencies

- **HiC-Pro** - Hi-C read processing
- **Bowtie2** - Read alignment
- **HiCExplorer** - Matrix normalization and visualization
- **Python 3.9+** with required libraries
- **SAMtools** and **BEDtools**

### Environment Configuration

```bash
# Create and activate conda environment for Hi-C analysis
conda create -n hic_analysis python=3.9
conda activate hic_analysis

# Install HiCExplorer
conda install -c bioconda hicexplorer

# Install Bowtie2 for alignment
conda install -c bioconda bowtie2

# Install additional utilities
conda install -c bioconda samtools bedtools
```

> **Note**: HiC-Pro may need to be installed manually from source if not available via conda.

---

## 🚀 Execution Steps

### 1. Prepare Genome Files

#### 1.1 Create Genome Size File

```bash
# Generate genome index if not exists
samtools faidx Genome.fasta

# Extract chromosome sizes
cut -f1,2 Genome.fasta.fai | sort -k1,1 -k2,2n > Genome.Size.txt
```

#### 1.2 Create Restriction Fragment File

```bash
# Generate DpnII restriction sites
/home/GuanJianTao/Software/HiC-Pro/bin/utils/digest_genome.py \
  -r DpnII \
  -o CG88.DpnII.bed \
  CG88.Chr.fa
```

---

### 2. Configure HiC-Pro

<details>
<summary><b>Click to view config-hicpro.txt configuration file</b></summary>

```bash
#########################################################################
## Paths and Settings  - Do not edit !
#########################################################################

TMP_DIR = /home/GuanJianTao/PanGenome.dir/HiC/hicpro/Genetic.Map.dir/tmp.dir/
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = rawdata

#######################################################################
## SYSTEM AND SCHEDULER - Start Editing Here !!
#######################################################################
N_CPU = 80
SORT_RAM = 2000M
LOGFILE = hicpro.log

JOB_NAME =
JOB_MEM =
JOB_WALLTIME =
JOB_QUEUE =
JOB_MAIL =

#########################################################################
## Data
#########################################################################

PAIR1_EXT = _R1
PAIR2_EXT = _R2

#######################################################################
## Alignment options
#######################################################################

MIN_MAPQ = 10

BOWTIE2_IDX_PATH = /home/GuanJianTao/PanGenome.dir/HiC/hicpro/Genetic.Map.dir/Bowtie2.Index.dir/
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = CG88.final.new.corrected
GENOME_SIZE = /home/GuanJianTao/PanGenome.dir/HiC/hicpro/Genetic.Map.dir/Bowtie2.Index.dir/Genome.Size.txt

#######################################################################
## Allele specific analysis
#######################################################################

ALLELE_SPECIFIC_SNP =

#######################################################################
## Capture Hi-C analysis
#######################################################################

CAPTURE_TARGET =
REPORT_CAPTURE_REPORTER = 1

#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = /home/GuanJianTao/PanGenome.dir/HiC/hicpro/Genetic.Map.dir/Bowtie2.Index.dir/CG88.DpnII.bed
LIGATION_SITE = GATC
MIN_FRAG_SIZE =
MAX_FRAG_SIZE =
MIN_INSERT_SIZE =
MAX_INSERT_SIZE =

#######################################################################
## Hi-C processing
#######################################################################

MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 0
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################

BIN_SIZE = 10000 20000 50000 100000 150000 200000 250000
MATRIX_FORMAT = upper

#######################################################################
## Normalization
#######################################################################
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1
```

</details>

---

### 3. Run HiC-Pro

```bash
/home/GuanJianTao/Software/HiC-Pro/bin/HiC-Pro \
  -i /home/GuanJianTao/PanGenome.dir/HiC/hicpro/Fastq/ \
  -o hic_pro.CG88.Map.out \
  -c config-hicpro.txt
```

---

### 4. Post-processing with HiCExplorer

#### 4.1 Convert HiC-Pro Matrices to H5 Format

```bash
# Define bin sizes
BIN_SIZES="10000 20000 50000 100000 150000 200000 250000"

for i in $BIN_SIZES
do
  hicConvertFormat \
    -m ../hicpro/Genetic.Map.dir/hic_pro.CG88.Map.out/hic_results/matrix/CG88/raw/${i}/CG88_${i}.matrix \
    --bedFileHicpro ../hicpro/Genetic.Map.dir/hic_pro.CG88.Map.out/hic_results/matrix/CG88/raw/${i}/CG88_${i}_abs.bed \
    --inputFormat hicpro \
    --outputFormat h5 \
    -o CG88.matrix.${i}.h5
done
```

#### 4.2 Apply ICE Correction

```bash
for i in $BIN_SIZES
do
  hicCorrectMatrix correct \
    --correctionMethod ICE \
    --matrix CG88.matrix.${i}.h5 \
    --filterThreshold -3.6 2 \
    -o CG88.corrected_matrix.${i}.h5
done
```

#### 4.3 Generate Diagnostic Plots (Optional)

```bash
hicCorrectMatrix diagnostic_plot \
  -m matrix.40000.h5 \
  -o matrix.40000.hic_corrected.png
```

#### 4.4 Plot Contact Matrices

```bash
for i in $BIN_SIZES
do
  hicPlotMatrix \
    -m CG88.corrected_matrix.${i}.h5 \
    -o CG88.matrix.${i}.hic_plot.whole.genome.pdf \
    --log1p \
    --dpi 300 \
    --colorMap YlOrRd \
    --chromosomeOrder chr1 chr2 chr3 chr4 chr5 chr6 chr7
done
```

---

## 📂 Output File Description

### HiC-Pro Outputs

| Directory | Description |
|-----------|-------------|
| `hic_results/` | Main results directory |
| `matrix/CG88/raw/` | Raw contact matrices at different resolutions |
| `matrix/CG88/iced/` | ICE-normalized matrices |
| `data/` | Processed alignment files |
| `pic/` | Quality control plots |

### HiCExplorer Outputs

| File Pattern | Description |
|--------------|-------------|
| `CG88.matrix.*.h5` | HiC-Pro matrices converted to H5 format |
| `CG88.corrected_matrix.*.h5` | ICE-corrected matrices |
| `CG88.matrix.*.hic_plot.whole.genome.pdf` | Contact matrix heatmaps |

---

## 📚 References

- [HiC-Pro Documentation](https://github.com/nservant/HiC-Pro)
- [HiCExplorer Documentation](https://hicexplorer.readthedocs.io/)
- [Hi-C Data Analysis Best Practices](https://www.nature.com/articles/nmeth.2148)

---

## 📧 Contact

For questions or issues, please open an issue on GitHub or contact the maintainer.

---

**Last Updated**: January 2026
