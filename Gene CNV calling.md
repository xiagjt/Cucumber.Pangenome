This code represents a **Pangenome Construction Pipeline** based on protein-to-genome alignments. It utilizes `miniprot` for high-speed mapping and `pangene` to build a gene-based pangenome graph.

---

# Pangenome Analysis Pipeline (miniprot & pangene)

This workflow describes how to construct a pangenome graph and generate a gene presence/absence (PAV) matrix using protein sequences mapped across multiple genome assemblies.

## 1. Prerequisites

The pipeline requires `miniprot`, `pangene`, and the `k8` Javascript runtime.

```bash
# Recommended: Install via Conda where available
conda create -n pangene_env -c bioconda miniprot pangene -y

# Note: k8 is often bundled with pangene or can be downloaded separately
# from https://github.com/attrth/k8

```

---

## 2. Protein-to-Genome Alignment

Map the representative protein sequences (pep) to each target genome assembly. This step must be repeated for all samples in your pangenome.

```bash
###### Map proteins to genome using miniprot
# --outs=0.97: Output alignments with score >= 97% of the maximum
# --no-cs: Do not output the cs tag
# -I: Indexing parameters
# -u: Un-spliced alignment (or specific splicing model)
# -t 16: Use 16 threads

miniprot --outs=0.97 --no-cs -Iut16 \
    CG110.final.new.fasta \
    CLv4.Gene.pep.Reprensentative.fa > CG110.paf

```

---

## 3. Pangenome Graph Construction

Once you have generated `.paf` files for all samples, combine them to build the pangenome graph.

```bash
###### Build the pangenome graph (GFA)
# This command takes all PAF files generated in the previous step
pangene *.paf > graph.gfa

```

---

## 4. Downstream Analysis (PAV and Variation)

Using the `k8` runtime and `pangene.js` script, you can extract meaningful biological data from the GFA graph.

### 4.1 Variation/Bubble Calling

Identify structural variations and PAVs (Presence/Absence Variations) within the graph.

```bash
###### Call "bubbles" (variants) from the graph
k8 pangene.js call graph.gfa > bubble.txt

```

### 4.2 Gene Presence/Absence Matrix

Generate an `.Rtab` file, which is a binary matrix (0/1) used for downstream comparative genomics or GWAS/GS analysis.

```bash
###### Generate Gene Presence/Absence (PAV) Matrix
k8 pangene.js gfa2matrix graph.gfa > gene_presence_absence.Rtab

```

---

This script provides a post-analysis step for your pangenome pipeline. After generating the **Gene Presence/Absence Matrix (`.Rtab`)** using `pangene`, you can use this Python script to classify genes into biological categories based on their copy number variation across the population.

---

## 5. Variant Classification (PAV & CNV)

This step categorizes genes from the pangenome matrix into four types:

* **PAV (Presence/Absence Variation):** Genes that are missing in at least one genome but present in others.
* **CNV (Copy Number Variation):** Genes that have more than one copy in at least one genome.
* **PAV & CNV:** Genes that exhibit both missingness and multi-copy events.
* **None (Core Single-copy):** Genes present as exactly one copy in all analyzed genomes.

### 5.1 Classification Script

<details>
<summary>Click to expand: classify_cnv_pav.py</summary>

```python
#!/usr/bin/env python3
"""
CNV and PAV Classification Script
Identifies Gene Copy Number Variation (CNV) and Presence/Absence Variation (PAV)
"""

import sys
import argparse
import pandas as pd
from pathlib import Path

def classify_variant(copy_numbers):
    """
    Classify variant type based on copy numbers.

    Parameters:
    -----------
    copy_numbers : array-like
        Gene copy numbers across different genomes.

    Returns:
    --------
    str : Variant type (PAV, CNV, PAV & CNV, or None)
    """
    try:
        copy_numbers = pd.to_numeric(copy_numbers, errors='coerce')
    except:
        return "Unknown"

    valid_numbers = copy_numbers.dropna()

    if len(valid_numbers) == 0:
        return "Unknown"

    # Check for presence of 0 (absence)
    has_zero = (valid_numbers == 0).any()
    has_non_zero = (valid_numbers > 0).any()

    # Check for copy number > 1
    has_cnv = (valid_numbers > 1).any()

    # Classification logic
    is_pav = has_zero and has_non_zero
    is_cnv = has_cnv

    if is_pav and is_cnv:
        return "PAV & CNV"
    elif is_pav:
        return "PAV"
    elif is_cnv:
        return "CNV"
    else:
        return "None"

def analyze_copy_number_matrix(input_file, output_file, genome_start_col=1,
                                genome_end_col=127, sep='\t'):
    """
    Analyze copy number matrix and classify CNV/PAV.
    """

    print(f"[INFO] Reading input file: {input_file}")

    try:
        df = pd.read_csv(input_file, sep=sep, low_memory=False)
    except Exception as e:
        print(f"[ERROR] Failed to read file: {e}")
        sys.exit(1)

    print(f"[INFO] Data dimensions: {df.shape[0]} rows x {df.shape[1]} columns")

    if genome_end_col > len(df.columns):
        print(f"[WARNING] Specified end column {genome_end_col} exceeds actual columns {len(df.columns)}")
        genome_end_col = len(df.columns)
        print(f"[WARNING] Adjusted to: {genome_end_col}")

    print(f"[INFO] Genome columns range: Index {genome_start_col} to {genome_end_col}")

    copy_number_cols = df.columns[genome_start_col:genome_end_col]

    print("[INFO] Starting classification...")

    variant_types = []
    for idx, row in df.iterrows():
        copy_numbers = row[copy_number_cols]
        variant_type = classify_variant(copy_numbers)
        variant_types.append(variant_type)

        if (idx + 1) % 1000 == 0:
            print(f"[PROGRESS] Processed {idx + 1} / {len(df)} genes...")

    df.insert(1, 'Variant_Type', variant_types)

    # Summary Statistics
    print("\n" + "="*60)
    print("Variant Type Statistics:")
    print("="*60)

    variant_counts = pd.Series(variant_types).value_counts()
    total_genes = len(variant_types)

    for variant_type, count in variant_counts.items():
        percentage = (count / total_genes) * 100
        print(f"{variant_type:20s}: {count:8d} ({percentage:6.2f}%)")

    print(f"{'Total':20s}: {total_genes:8d} (100.00%)")
    print("="*60 + "\n")

    # Save Results
    print(f"[INFO] Saving results to: {output_file}")
    try:
        df.to_csv(output_file, sep=sep, index=False)
        print("[SUCCESS] Results saved successfully!")
    except Exception as e:
        print(f"[ERROR] Failed to save file: {e}")
        sys.exit(1)

    # Save Stats
    stats_file = output_file.rsplit('.', 1)[0] + '_statistics.txt'
    with open(stats_file, 'w') as f:
        f.write("CNV and PAV Classification Statistics\n")
        f.write("="*60 + "\n")
        f.write(f"Total genes: {total_genes}\n")
        for variant_type, count in variant_counts.items():
            percentage = (count / total_genes) * 100
            f.write(f"{variant_type:20s}: {count:8d} ({percentage:6.2f}%)\n")

    return df, variant_counts

def main():
    parser = argparse.ArgumentParser(
        description='Classify Gene PAV and CNV from copy number matrix',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-i', '--input', required=True, help='Input matrix path')
    parser.add_argument('-o', '--output', required=True, help='Output path')
    parser.add_argument('-s', '--start-col', type=int, default=1, help='Start index of genome columns (0-based, default=1)')
    parser.add_argument('-e', '--end-col', type=int, default=127, help='End index of genome columns (default=127)')
    parser.add_argument('--sep', type=str, default='\t', help='Separator (default: tab)')

    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"[ERROR] Input file not found: {args.input}")
        sys.exit(1)

    analyze_copy_number_matrix(
        input_file=args.input,
        output_file=args.output,
        genome_start_col=args.start_col,
        genome_end_col=args.end_col,
        sep=args.sep
    )

if __name__ == '__main__':
    main()

```

</details>

### 5.2 Usage Instructions

Ensure you have the `pandas` library installed in your environment before running.

```bash
# Install dependency
pip install pandas

# Run classification
# Assuming your gene_presence_absence.Rtab has genome data from column 2 to 127
python classify_cnv_pav.py \
    -i gene_presence_absence.Rtab \
    -o pangenome_classified.txt \
    -s 1 \
    -e 127

```

**Key Outputs:**

* `pangenome_classified.txt`: The original matrix with an added `Variant_Type` column.
* `pangenome_classified_statistics.txt`: A summary report of the PAV/CNV distribution.

---

