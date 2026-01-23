This section details the **Selection Analysis** used to identify selective sweeps and genomic regions under pressure during cucumber evolution or domestication. We employ two complementary statistics: **** (Fixation Index) and **XP-CLR** (Cross-Population Composite Likelihood Ratio).

---

# Selection Analysis (Selective Sweep Detection)

## 1. Environment Setup

We utilize **VCFtools** for  calculations and **Bedtools** for genomic window merging.

```bash
# Create environment for selection analysis
conda create -n selection_env -c bioconda vcftools bedtools -y

# Activate the environment
conda activate selection_env

```

---

## 2.  (Fixation Index) Analysis

 is used to measure genetic differentiation between populations. High  values indicate regions potentially under divergent selection.

### 2.1 Sliding Window Calculation

We calculate  using a  window with a  step size between specified populations (e.g., Indian vs. XIS).

```bash
###### Calculate Fst between populations
vcftools --vcf ./merge.InDel.SNP.SV.Anno.422.vcf \
    --weir-fst-pop ./Indian.Samples.list \
    --weir-fst-pop ./XIS.Samples.list \
    --fst-window-size 20000 \
    --fst-window-step 2000 \
    --out ./Fst.Indian.XIS.out

```

### 2.2 Filtering and Extracting Top Hits

We filter windows with low variant density (at least 5 variants) and extract the top 5% of differentiated regions.

```bash
# 1. Filter for variant density and sort by Fst value
for f in *.fst
do
    awk '{if($4>=5){print $0}}' ${f} | sort -k5 -g -r > ${f}.sort.tsv 
done 

# 2. Extract top 5% windows (Replace XXXX with 5% of total window count)
head -n XXXX ${f}.sort.tsv > ${f}.sort.top5.tsv

# 3. Merge overlapping windows using Bedtools
for f in *.fst
do
    sort -k1,1 -k2,2n ${f}.sort.top5.tsv > ${f}.sort.top5.sort.tsv
    bedtools merge -c 5 -o mean -i ${f}.sort.top5.sort.tsv > ${f}.sort.top5.merge.tsv
done

```

---

## 3. XP-CLR Analysis

XP-CLR is a multi-locus test that uses the allele frequency spectrum to detect selective sweeps by comparing two populations.

### 3.1 Processing Multi-Population Comparisons

The following loops process XP-CLR results across multiple population pairs (e.g., East Asian vs. European, Indian vs. XIS).

```bash
###### Concatenate and filter XP-CLR results
for i in EA.vc.EU EA.vc.Indian EA.vc.XIS EU.vc.Indian Indian.vc.XIS EU.vc.XIS
do
    # Merge chromosome results and filter for score >= 5
    cat ${i}.chr*.out | egrep -v "^id" | \
    awk '{if($10>=5){print $0}}' | \
    sort -k12 -g -r > ${i}.xpclr.genome.sort.out
    
    # Extract top 5% hits
    head -n XXXX ${i}.xpclr.genome.sort.out > ${i}.xpclr.genome.sort.top5.out
done

```

### 3.2 Genomic Window Merging

To identify "Selection Signals" instead of individual windows, we merge adjacent high-scoring regions.

```bash
###### Merge XP-CLR windows using Bedtools
for i in EA.vc.EU EA.vc.Indian EA.vc.XIS EU.vc.Indian Indian.vc.XIS EU.vc.XIS
do
    # Prepare BED-like format (Chr, Start, End, Score)
    sort -k2,2 -k3,3n ${i}.xpclr.genome.sort.top5.out | \
    cut -f 2-4,12 > ${i}.xpclr.genome.sort.top5.sort.out
    
    # Merge overlapping windows and calculate mean XP-CLR score
    bedtools merge -c 4 -o mean -i ${i}.xpclr.genome.sort.top5.sort.out > ${i}.xpclr.genome.sort.top5.merge.out
done

```

---
