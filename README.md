# Cucumber Pangenome Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains custom scripts and computational pipelines for pangenome analysis, structural variant detection, and multi-omics integration in cucumber genomics research.

## Citation

If you use this code in your research, please cite:
```
Under review
```

## Table of Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Analysis Scripts](#analysis-scripts)
- [Usage Examples](#usage-examples)
- [Input and Output Files](#input-and-output-files)
- [License](#license)
- [Contact](#contact)

## Overview

This pipeline integrates multiple bioinformatics tools for comprehensive cucumber pangenome analysis, including:

- **Genome assembly and quality assessment** from PacBio HiFi long reads
- **Structural variant (SV) detection** and graph-based pangenome construction
- **SNP/InDel calling and annotation** from whole-genome resequencing data
- **Gene and repeat annotation** with functional classification
- **Syntenic pan-genome gene analysis** including core/dispensable gene identification
- **Gene-level PAV and CNV detection** across populations
- **Population genetics** including phylogenetics, FST, and selection scans
- **GWAS and eQTL analysis** for trait-associated variant identification
- **TWAS and colocalization** integrating genomics with transcriptomics
- **NLR resistance gene annotation** and classification
- **Haplotype diversity analysis** and group-specific variant identification

## System Requirements

### Software Dependencies

All software used in this pipeline are open-source tools with OSI-approved or permissive licenses. Major tools include:

#### Assembly & Alignment
- Hifiasm v0.19.5
- RagTag v2.1.0
- Minimap2 v2.17-r941 / v2.28-r1209
- BWA v0.7.17-r1198
- pbmm2 v1.13.1
- MUMmer v4.0

#### Variant Calling
- GATK v4.1.5.0
- Sniffles v2.3.2
- SVIM v2.0.0
- CuteSV v2.1.1
- PBSV v2.9.0
- Syri
- SURVIVOR v1.0.6
- vg toolkit
- Paragraph

#### Quality Assessment
- BUSCO v5.2.2
- Merqury v1.1
- fastp

#### Gene Annotation
- RepeatModeler v2.0.3
- RepeatMasker v4.1.2-p1
- Hisat2 v2.2.1
- STAR v2.7.10a
- StringTie v2.2.1
- EVidenceModeler v1.1.1
- BRAKER2
- GenomeThreader
- InterProScan
- DIAMOND v2.0.15

#### Statistical Analysis & Association Studies
- VCFtools v0.1.17
- PLINK
- Beagle v21Jan17.6cc
- GEMMA v0.98.5
- EMMAX v2012-02-10 / beta-07Mar2010
- COLOC v5.2.3
- R/qtl
- LepMap3
- SNPbinner

#### Phylogenetics & Population Structure
- IQ-TREE v2.3.6
- GCTA v1.94.1
- ADMIXTURE v1.3.0

#### Pan-genome Analysis
- mSynOrths v0.1 (https://gitee.com/zhanglingkui/msynorths)
- PAML (CODEML)
- ClustalW2
- PAL2NAL
- miniprot
- pangene
- CD-HIT v4.8.1
- GMAP

#### Chromatin & RNA Analysis
- HiC-Pro
- HiCExplorer
- DESeq2 v1.30.1 (R)
- clusterProfiler v3.18.1 (R)

#### Utilities
- Picard v1.118
- BEDTools
- TelomereSearch (https://github.com/jamiemcg/TelomereSearch)
- NLR-annotator

#### Programming Languages
- Python >= 3.8
- R >= 4.0
- Bash shell

### Operating System
- **Tested on:** Linux (Ubuntu 20.04 LTS, CentOS 7)
- **Compatible with:** Most Unix-like systems
- **Note:** MacOS should work but not fully tested

### Hardware Requirements

#### Minimum Requirements
- **RAM:** 32 GB
- **CPU:** 8 cores
- **Storage:** 500 GB free disk space

#### Recommended for Full Dataset Analysis
- **RAM:** 128 GB or more
- **CPU:** 32+ cores
- **Storage:** 2 TB+ free disk space
- **Note:** Most analyses (especially genome assembly, SV calling, and GWAS) are computationally intensive and benefit from high-performance computing (HPC) environments

## Installation

### Software Installation

All bioinformatics tools used in this pipeline can be installed via:

1. **Conda/Mamba (Recommended):** Most tools are available through Bioconda channel
```bash
conda install -c bioconda <tool_name>
```

2. **Manual Installation:** Follow the installation instructions provided by each tool's author in their respective repositories/documentation

### Clone This Repository
```bash
git clone https://github.com/xiagjt/Cucumber.Pangenome.git
cd Cucumber.Pangenome
```

## Analysis Scripts

This repository contains scripts for the following analyses:

### 1. Genome Assembly and Anchoring
**Script:** `Genome assembly and anchoring`

**Purpose:** Generate chromosome-level genome assemblies from PacBio HiFi reads

**Key Steps:**
- Genome assembly using Hifiasm
- Chromosome anchoring with RagTag to reference genome
- Quality assessment using BUSCO, Merqury, and read alignment validation

---

### 2. Gene Annotation
**Script:** `Gene annotation`

**Purpose:** Predict protein-coding genes and assign functional annotations

**Key Steps:**
- Gene prediction integration using EVidenceModeler (combining BRAKER2, GenomeThreader, StringTie)
- Functional annotation with InterProScan and DIAMOND (NR, GO, KEGG databases)

---

### 3. Repeat Annotation
**Script:** `Repeat annotation`

**Purpose:** Identify and classify repetitive elements in genomes

**Key Steps:**
- De novo repeat library construction with RepeatModeler
- Repeat annotation with RepeatMasker using custom and Repbase libraries

---

### 4. Centromere and Telomere Detection
**Script:** `Centromere and Telomere detection`

**Purpose:** Identify centromeric and telomeric regions

**Key Steps:**
- Centromere identification based on 177-bp satellite repeat enrichment
- Telomere detection using TelomereSearch

---

### 5. SNP Calling and Annotation
**Script:** `SNP calling and annotation`

**Purpose:** Call and annotate single nucleotide polymorphisms and small InDels

**Key Steps:**
- Read alignment with BWA-MEM
- Variant calling with GATK HaplotypeCaller
- Variant filtering and quality control
- Functional annotation of variants

---

### 6. SV Calling and Graph-Pangenome Construction
**Script:** `SV Calling and graph-Pangenome construction`

**Purpose:** Detect structural variants and construct graph-based pangenome

**Key Steps:**
- Multi-caller SV detection (Sniffles, SVIM, CuteSV, PBSV, Syri)
- SV merging across callers using SURVIVOR
- Graph-based pangenome construction with vg toolkit
- SV genotyping across populations using short-read data

---

### 7. FST and XP-CLR Analysis
**Script:** `Fst and XP-CLR`

**Purpose:** Identify genomic regions under selection

**Key Steps:**
- FST calculation between populations using VCFtools
- XP-CLR analysis for detecting selective sweeps
- Identification of candidate selective regions (top 5% empirical distribution)

---

### 8. Gene-based Pangenome Analysis
**Script:** `Syntenic pangenome`

**Purpose:** Construct syntenic pangenome and classify gene families as well as haplotype analysis

**Key Steps:**
- Syntenic gene family clustering using mSynOrths
- Gene family classification (core/softcore/dispensable/private)
- Ka/Ks calculation with PAML CODEML
- Haplotype analysis and group-specific frequency testing

---

### 9. Specific and Shared Gene Families
**Script:** `Specific and shared gene families`

**Purpose:** Identify lineage-specific and shared genes across groups

**Key Steps:**
- Detection of group-specific haplotypes
- Chi-square testing for frequency differences (FDR < 0.001)
- Identification of genes unique to pangenome vs. reference

---

### 10. Gene CNV Calling
**Script:** `Gene CNV calling`

**Purpose:** Detect gene-level copy number variations

**Key Steps:**
- Protein alignment to genomes using miniprot
- CNV identification with pangene tool
- Population differentiation analysis using VST statistic
- Redundancy removal with CD-HIT

---


### 11. NLR Annotation
**Script:** `NLR annotation`

**Purpose:** Identify and classify nucleotide-binding leucine-rich repeat (NLR) resistance genes

**Key Steps:**
- Initial NLR identification using InterProScan
- Novel NLR prediction with BRAKER2 and NLR-annotator
- NLR classification (TNLs, CNLs, RNLs, NLs) based on domain composition

---

### 12. GWAS and eQTL Analysis
**Script:** `GWAS and eQTL`

**Purpose:** Perform genome-wide association studies and expression QTL mapping

**Key Steps:**
- SNP/SV-based GWAS with GEMMA (mixed linear model)
- eQTL mapping with EMMAX (SV-expression associations)
- Cis/trans-eQTL classification
- eQTL hotspot detection

---

### 13. TWAS and GWAS-eQTL Colocalization
**Script:** `TWAS and GWAS-eQTL Colocalization`

**Purpose:** Integrate transcriptome with genome-wide association results

**Key Steps:**
- Transcriptome-wide association study (TWAS) using EMMAX
- Colocalization analysis with COLOC
- Identification of causal variants (H4 posterior probability > 0.50)

### 14. Hi-C heatmap plot of CG88 for validating the large inversions
**Script:** `Hi-C heatmap`

---


## Usage Examples

### Example 1: Genome Assembly
```bash
# Navigate to the genome assembly script
cd "Genome assembly and anchoring/"
# Follow the instructions in the script for Hifiasm assembly and RagTag anchoring
```

### Example 2: SNP Calling
```bash
# Navigate to SNP calling script
cd "SNP calling and annotation/"
# Run the SNP calling pipeline according to script instructions
```

### Example 3: SV Detection and Pangenome Construction
```bash
# Navigate to SV calling script
cd "SV Calling and graph-Pangenome construction/"
# Run multi-caller SV detection and graph pangenome construction
```

### Example 4: Gene Family Analysis
```bash
# Navigate to gene-based pangenome script
cd "Gene-based pangenome/"
# Run syntenic gene family clustering and Ka/Ks analysis
```

### Example 5: GWAS and eQTL
```bash
# Navigate to GWAS and eQTL script
cd "GWAS and eQTL/"
# Run GWAS and eQTL mapping analyses
```

### Example 6: NLR Gene Annotation
```bash
# Navigate to NLR annotation script
cd "NLR annotation/"
# Run NLR gene annotation pipeline
```

## Input and Output Files

### Common Input File Formats

| File Type | Format | Description |
|-----------|--------|-------------|
| Genome sequences | FASTA (`.fa`, `.fasta`) | Reference genome or assembled contigs |
| Sequencing reads | FASTQ (`.fq`, `.fastq`, `.fq.gz`) | Raw sequencing data |
| Aligned reads | BAM (`.bam`) | Aligned and sorted reads with index |
| Variants | VCF (`.vcf`, `.vcf.gz`) | Variant calls (SNPs, InDels, SVs) |
| Gene annotations | GFF3/GTF | Gene structure annotations |
| Phenotype data | Tab-delimited TXT | Sample ID and trait values |
| Expression data | Tab-delimited TXT | Gene expression matrix (TPM/FPKM) |

### Key Output Files

| Analysis Script | Output File | Description |
|-----------------|-------------|-------------|
| Genome assembly | `*.assembly.fasta` | Chromosome-level genome assembly |
| Assembly QC | `busco_summary.txt` | Gene completeness statistics |
| Gene annotation | `*.annotation.gff3` | Predicted gene structures |
| Gene annotation | `functional_annotation.txt` | GO, KEGG assignments |
| SNP calling | `*.snps.vcf.gz` | Filtered SNPs and InDels |
| SV calling | `*.svs.vcf.gz` | Structural variant calls |
| Pangenome | `pangenome_graph.vg` | Graph-based pangenome |
| Gene families | `gene_families.txt` | Core/softcore/dispensable/private genes |
| Gene families | `ka_ks_results.txt` | Selection pressure analysis |
| Gene CNV | `gene_cnv_matrix.txt` | Gene copy number matrix |
| FST analysis | `fst_results.txt` | FST values per window |
| XP-CLR analysis | `xpclr_results.txt` | XP-CLR scores per window |
| GWAS | `gwas_results.txt` | Association P-values |
| eQTL | `eqtl_results.txt` | Expression-associated variants |
| TWAS | `twas_results.txt` | Expression-trait associations |
| Colocalization | `coloc_results.txt` | H4 posterior probabilities |
| NLR annotation | `nlr_genes.gff3` | Annotated NLR resistance genes |

## Repository Structure
```
Cucumber.Pangenome/
├── Centromere and Telomere detection
├── Fst and XP-CLR
├── GWAS and eQTL
├── Gene CNV calling
├── Gene annotation
├── Gene-based pangenome
├── Genome assembly and anchoring
├── Haplotype diversity
├── Haplotype frequency comparison
├── Homologous group determination
├── NLR annotation
├── Repeat annotation
├── SNP calling and annotation
├── SV Calling and graph-Pangenome construction
├── Specific and shared gene families
├── TWAS and GWAS-eQTL Colocalization
├── LICENSE
└── README.md
```

## Notes

- Most analyses require substantial computational resources; HPC is recommended for full-scale studies
- Intermediate files can be large; ensure sufficient disk space
- Scripts assume input files follow standard bioinformatics formats
- Modify parameters in scripts according to your system capabilities and data characteristics
- Refer to comments within each script for detailed parameter descriptions and usage instructions

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this pipeline in your research, please cite:
```
under review
```

## Contact

For questions and feedback:
- **GitHub Issues:** https://github.com/xiagjt/Cucumber.Pangenome/issues
- **Corresponding Author:** Jiantao Guan, guanjiantao@caas.cn

## Acknowledgments

We thank the developers of all open-source bioinformatics tools used in this pipeline. This work was supported by the National Natural Science Foundation of China and other funding agencies.

---

**Last Updated:** January 2026
