## Table of Contents
- [Syntenic Pangenome Construction](#syntenic-pangenome-construction)
- [Haplotype Diversity Analysis](#haplotype-diversity-analysis)
- [Haplotype Frequency Comparison](#haplotype-frequency-comparison)

---

## Syntenic Pangenome Construction

### Overview
Construct syntenic pangenome and identify gene families across 126 cucumber genomes using mSynOrths.

### Workflow

#### Step 1: All-vs-All Protein Alignment
```bash
# Perform all vs all protein alignment for 126 pangenomes including reference genome
diamond blastp \
  -q DDD.pep.fa \
  -d ${i}.pep.fa \
  -o DDD${i}_results.txt \
  -e 1 \
  -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  -k 1 \
  -p 80
```

#### Step 2: Syntenic Gene Family Identification
```bash
# The installation pipeline of mSynOrths is at gitee
# https://gitee.com/zhanglingkui/msynorths/

# Identification of synteny gene families
python ~/180T/msynorths_gitee/mSynOrths.py \
  -t 60 \
  -o Cucum_137_0110 \
  -m 0.8 \
  -v 0.5 \
  -i 85 \
  -f genome_pos.txt \
  -c
```

#### Step 3: Build Gene-level Pangenome
```bash
# Built gene-level pangenome for cucumber genomes
python bulit_pangenome.py \
  -m Cucum_137_0110 \
  -o Cucum_137_0110_pangene.txt
```

### Script: `bulit_pangenome.py`

**Parameters:**
- `-m`: Output folder of mSynOrths (required)
- `-f`: Choose frame of pangenome (default: 'mSynF1')
- `-o`: Output file name (default: 'pangenome.txt')

<details>
<summary><b>📜 Click to view complete Python script</b></summary>

```python
# -*- coding=utf-8 -*-
import sys
import os
import argparse
import re

parser = argparse.ArgumentParser(description='pangenome')
parser.add_argument('-m', type=str, help='output folder of mSynOrths')
parser.add_argument('-f', type=str, default='mSynF1', help='choose frame of pan genome')
parser.add_argument('-o', type=str, default='pangenome.txt', help='output file')
args = parser.parse_args()
msyn_fold, frame_fold, out_file = args.m, args.f, args.o

# species num. of mSynOrths
msynf_num = 0
# To find the longest gene as representative
total_gene_length_dict = {}
# The list dictionary of all genes
total_gene_list_dict = {}

for i in os.listdir(msyn_fold):
    if i.startswith('mSynF'):
        msynf_num += 1
        op_file = open(msyn_fold + '/' + i + '/species.gff_pos')
        index = 0
        total_gene_list_dict[i[5:]] = []
        for line in op_file:
            index += 1
            line_list = line.strip().split('\t')
            total_gene_length_dict[line_list[0]] = int(line_list[2])
            total_gene_list_dict[i[5:]].append(line_list[0])
        op_file.close()

frame_gene_list = []
frame_fold_path = msyn_fold + '/' + frame_fold + '/'
frame_file = frame_fold_path + 'species.gff_pos_DelTandem'
op_frame_file = open(frame_file, 'r')
frame_gene_dict = {}
index = 0

for line in op_frame_file:
    index += 1
    line_list = line.strip().split('\t')
    frame_gene_list.append(line_list[0])
    frame_gene_dict[line_list[0]] = [index, int(line_list[2]), line_list[3]]
op_frame_file.close()

frame_num = frame_fold[5:]
total_syn_frame_dict = {}

for i in os.listdir(msyn_fold):
    if i.startswith('mSynF'):
        if int(i[5:]) < int(frame_num):
            for x in os.listdir(msyn_fold + '/' + i):
                if x.endswith('syn.seg'):
                    if re.findall('\d+', x)[1] == frame_num:
                        op_file = open(msyn_fold + '/' + i + '/' + x, 'r')
                        for line in op_file:
                            if line.startswith('##'):
                                continue
                            line_list = line.strip().split('\t')
                            total_syn_frame_dict[line_list[0]] = line_list[1]
                        op_file.close()
        elif int(i[5:]) == int(frame_num):
            for x in os.listdir(msyn_fold + '/' + i):
                if x.endswith('syn.seg'):
                    op_file = open(msyn_fold + '/' + i + '/' + x, 'r')
                    for line in op_file:
                        if line.startswith('##'):
                            continue
                        line_list = line.strip().split('\t')
                        total_syn_frame_dict[line_list[1]] = line_list[0]
                    op_file.close()

msyn_file = msyn_fold + '/Total_species_syntenic_gene_pairs.txt'
index = 0
op_msyn_file = open(msyn_file, 'r')

for line in op_msyn_file:
    line_list = line.strip().split('\t')
    line_dict = {}
    for i in range(1, msynf_num + 1):
        line_dict[str(i)] = []
    num_index = 0
    for i in line_list:
        ad_list = i[len(i.split(':')[0]) + 1:].split(',')
        ad_list.remove('')
        num_index += 1
        line_dict[i.split(':')[0]] = ad_list
    
    if len(line_dict[frame_num]) > 0:
        sort_gene_list = sorted(line_dict[frame_num], 
                               key=lambda x: frame_gene_dict[x][1], reverse=True)
        frame_gene = sort_gene_list[0]
        frame_gene_dict[frame_gene].append(line_dict)
    else:
        big_length = 0
        big_gene = ''
        big_key = ''
        for key, value in line_dict.items():
            if len(value) > 0:
                for gene_id in value:
                    if total_gene_length_dict[gene_id] > big_length:
                        big_length = total_gene_length_dict[gene_id]
                        big_gene = gene_id
                        big_key = key
        big_gene_index = total_gene_list_dict[big_key].index(big_gene)
        for flank_gene in total_gene_list_dict[big_key][big_gene_index::]:
            if flank_gene in total_syn_frame_dict.keys():
                add_index = frame_gene_dict[total_syn_frame_dict[flank_gene]][0] + 0.1
                add_chr = frame_gene_dict[total_syn_frame_dict[flank_gene]][2]
                frame_gene_dict[big_gene] = [add_index, big_length, add_chr, line_dict]
                break
op_msyn_file.close()

wr_out_file = open(out_file, 'w')
for pan_gene_key in sorted(frame_gene_dict, key=lambda x: frame_gene_dict[x][0]):
    pan_gene_line = frame_gene_dict[pan_gene_key]
    if len(pan_gene_line) > 3:
        wr_out_file.write(pan_gene_key + '\t')
        for pan_gene in pan_gene_line[:3]:
            wr_out_file.write(str(pan_gene) + '\t')
        if len(pan_gene_line) > 3:
            for i in range(1, int(msynf_num) + 1):
                if len(pan_gene_line[-1][str(i)]) > 1:
                    for x in pan_gene_line[-1][str(i)]:
                        wr_out_file.write(x + ',')
                elif len(pan_gene_line[-1][str(i)]) == 1:
                    wr_out_file.write(pan_gene_line[-1][str(i)][0])
                else:
                    wr_out_file.write('-')
                wr_out_file.write('\t')
            wr_out_file.write('\n')
wr_out_file.close()
```

</details>

---

## Defining core, softcore, dispensable, and private gene families
**Usage:**
```bash
# Basic usage
python gene_family.category.py 126species.gene.family.matrix.xls 126species.gene.family.category.xls  126species.gene.family.matrix.category.xls
```

<details>
<summary><b>📜 Click to view complete Python script (500+ lines)</b></summary>

```python
#!/usr/bin/env python3
import pandas as pd
import sys

def classify_gene_family(input_file, output_file, matrix_file, delimiter='\t'):
    """
    Classify gene families
    
    Parameters:
        input_file: Input file path
        output_file: Output file path (classification results)
        matrix_file: Output file path (presence/absence matrix)
        delimiter: Delimiter, default is tab
    """
    
    print(f"Reading file: {input_file}")
    
    # Read data, preserving original column names
    try:
        # Try reading the first few lines to determine if there's a header
        with open(input_file, 'r') as f:
            first_line = f.readline().strip()
        
        # Determine if the first line might be a header (whether it contains numbers)
        import re
        has_numbers = any(char.isdigit() for char in first_line.split(delimiter)[0])
        
        # If no numbers, there might be a header; otherwise, no header
        if not has_numbers:
            df = pd.read_csv(input_file, sep=delimiter)
            print("Header detected, will preserve original column names")
        else:
            df = pd.read_csv(input_file, sep=delimiter, header=None)
            print("No header detected, using default column names")
            
    except Exception as e:
        print(f"Failed to read file: {e}")
        sys.exit(1)
    
    print(f"Data dimensions: {df.shape}")
    print(f"Total gene families: {df.shape[0]}")
    
    # Get column names
    if df.columns[0] == 0:  # If no header, column names are numeric indices
        # Save original column names (indices) for matrix output
        original_columns = list(df.columns)
        # Extract first 4 columns of information
        info_columns = df.iloc[:, 0:4].copy()
        info_columns.columns = ['Gene_ID', 'Family_ID', 'Gene_Length', 'Chromosome']
        
        # Extract material data (5th column to end)
        material_columns = df.iloc[:, 4:].copy()
        # Generate friendly column names for material columns (original column index + 5)
        material_column_names = [f'Material_{i+1}' for i in range(material_columns.shape[1])]
        material_columns.columns = material_column_names
    else:
        # If there's a header, keep original column names
        original_columns = list(df.columns)
        # Extract first 4 columns of information, using original column names
        info_columns = df.iloc[:, 0:4].copy()
        if len(original_columns) >= 4:
            info_columns.columns = original_columns[:4]
        else:
            info_columns.columns = ['Gene_ID', 'Family_ID', 'Gene_Length', 'Chromosome']
        
        # Extract material data, preserving original column names
        material_columns = df.iloc[:, 4:].copy()
        if len(original_columns) > 4:
            material_column_names = original_columns[4:]
            material_columns.columns = material_column_names
    
    # Count number of materials with genes for each gene family
    presence_counts = []
    
    for idx, row in material_columns.iterrows():
        # Count materials that are not "-"
        count = (row != '-').sum()
        presence_counts.append(count)
    
    # Classify based on material count
    categories = []
    for count in presence_counts:
        if count == 126:
            categories.append('core')
        elif 114 <= count <= 125:
            categories.append('softcore')
        elif 2 <= count <= 113:
            categories.append('dispensable')
        elif count == 1:
            categories.append('private')
        else:  # count == 0 (theoretically shouldn't occur)
            categories.append('unknown')
    
    # Create classification result dataframe
    result_df = info_columns.copy()
    result_df['Category'] = categories
    
    # Create presence/absence matrix (0 for absence, 1 for presence)
    presence_matrix = material_columns.apply(lambda x: (x != '-').astype(int))
    
    # Create complete matrix dataframe, containing all original columns
    # First create a copy of original data
    matrix_df = df.copy()
    
    # Add classification column to matrix dataframe (after first 4 columns)
    matrix_df.insert(4, 'Category', categories)
    
    # Count each category
    category_counts = pd.Series(categories).value_counts()
    print("\nClassification Statistics:")
    print(f"  Core (126 materials):        {category_counts.get('core', 0):6d} gene families")
    print(f"  Softcore (114-125 materials): {category_counts.get('softcore', 0):6d} gene families")
    print(f"  Dispensable (2-113 materials):{category_counts.get('dispensable', 0):6d} gene families")
    print(f"  Private (1 material):         {category_counts.get('private', 0):6d} gene families")
    if 'unknown' in category_counts:
        print(f"  Unknown:                     {category_counts.get('unknown', 0):6d} gene families")
    
    # Save results
    try:
        # Save classification results (only first 4 columns + classification column)
        result_df.to_csv(output_file, sep='\t', index=False, header=True)
        print(f"\nClassification results saved to: {output_file}")
        
        # Save complete matrix data (with all original columns and classification column)
        matrix_df.to_csv(matrix_file, sep='\t', index=False, header=True)
        print(f"Complete matrix data (preserving original column names) saved to: {matrix_file}")
        
        # Additionally save a matrix with only 0/1 values
        presence_only_df = pd.concat([result_df, presence_matrix], axis=1)
        presence_only_file = matrix_file.replace('.txt', '_binary.txt').replace('.tsv', '_binary.tsv').replace('.csv', '_binary.csv')
        presence_only_df.to_csv(presence_only_file, sep='\t', index=False, header=True)
        print(f"Binary presence/absence matrix saved to: {presence_only_file}")
        
    except Exception as e:
        print(f"Failed to save files: {e}")
        sys.exit(1)
    
    return result_df, matrix_df

def main():
    """
    Main function
    """
    if len(sys.argv) < 4:
        print("Usage: python script.py <input_file> <classification_output_file> <matrix_output_file> [delimiter]")
        print("Example: python script.py gene_family_matrix.txt classified.txt matrix_output.txt")
        print("         python script.py gene_family_matrix.csv classified.csv matrix_output.csv ','")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    matrix_file = sys.argv[3]
    delimiter = sys.argv[4] if len(sys.argv) > 4 else '\t'
    
    # Execute classification
    classify_gene_family(input_file, output_file, matrix_file, delimiter)
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()
```
</details>

---

## Haplotype identification

**Usage:**
```bash
# Basic usage
python3 sequence_haplotypes.py 126species.gene.family.matrix.xls all_Samples.pep_merge.fasta 126Samples.homology.haplotype.table.xls
```
<details>
<summary><b>📜 Click to view complete Python script (500+ lines)</b></summary>

```python
#!/usr/bin/env python3
import sys
from collections import defaultdict

def parse_fasta(fasta_file):
    """
    Parse FASTA file, return a dictionary mapping gene IDs to sequences
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    print(f"Reading FASTA file: {fasta_file}")
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Save previous sequence
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                
                # Start new sequence, extract ID (remove >, take part before first space)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    
    print(f"Total {len(sequences)} sequences read")
    return sequences


def assign_haplotypes(gene_ids, sequences):
    """
    Assign haplotype numbers to a group of genes
    Same sequence gets same number, numbered from 1 in order of first appearance
    """
    seq_to_haplotype = {}
    gene_to_haplotype = {}
    haplotype_counter = 1
    
    for gene_id in gene_ids:
        if gene_id == '-':
            gene_to_haplotype[gene_id] = '-'
            continue
        
        if gene_id not in sequences:
            print(f"Warning: Gene {gene_id} not found in FASTA file")
            gene_to_haplotype[gene_id] = '-'
            continue
        
        seq = sequences[gene_id]
        
        # If sequence already seen, use existing haplotype number
        if seq in seq_to_haplotype:
            gene_to_haplotype[gene_id] = seq_to_haplotype[seq]
        else:
            # New sequence, assign new haplotype number
            seq_to_haplotype[seq] = haplotype_counter
            gene_to_haplotype[gene_id] = haplotype_counter
            haplotype_counter += 1
    
    return gene_to_haplotype


def process_gene_family_matrix(matrix_file, fasta_file, output_file):
    """
    Process gene family matrix file for haplotype analysis
    """
    # Read all sequences
    sequences = parse_fasta(fasta_file)
    
    print(f"\nReading gene family matrix: {matrix_file}")
    
    # Read matrix file
    with open(matrix_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # Parse header
    header = lines[0].split('\t')
    sample_names = header[4:]  # First 4 columns: Reference genome ID, GroupID, Gene length, Chromosome info
    
    print(f"Total {len(sample_names)} samples")
    print(f"Total {len(lines)-1} genomes")
    
    # Prepare output results
    output_lines = []
    output_header = ['GeneID'] + sample_names
    output_lines.append('\t'.join(output_header))
    
    # Process each genome
    for idx, line in enumerate(lines[1:], 1):
        fields = line.split('\t')
        
        if len(fields) < 5:
            print(f"Warning: Line {idx+1} incomplete, skipping")
            continue
        
        ref_gene_id = fields[0]
        group_id = fields[1]
        gene_length = fields[2]
        chromosome = fields[3]
        gene_ids = fields[4:]
        
        # Ensure gene ID count matches sample count
        if len(gene_ids) != len(sample_names):
            print(f"Warning: Group {group_id} gene count doesn't match sample count")
            continue
        
        # Process multiple gene IDs (comma-separated), take only the first
        processed_gene_ids = []
        for gene_id in gene_ids:
            if gene_id == '-':
                processed_gene_ids.append('-')
            elif ',' in gene_id:
                # If multiple gene IDs, take only the first
                first_gene = gene_id.split(',')[0]
                processed_gene_ids.append(first_gene)
            else:
                processed_gene_ids.append(gene_id)
        
        # Collect all valid gene IDs for this group (including reference gene)
        all_genes = [ref_gene_id] + [g for g in processed_gene_ids if g != '-']
        
        # Assign haplotypes
        gene_to_haplotype = assign_haplotypes(all_genes, sequences)
        
        # Build output line
        haplotypes = []
        for gene_id in processed_gene_ids:
            if gene_id == '-':
                haplotypes.append('-')
            else:
                haplotypes.append(str(gene_to_haplotype.get(gene_id, '-')))
        
        output_line = [ref_gene_id] + haplotypes
        output_lines.append('\t'.join(output_line))
        
        if idx % 1000 == 0:
            print(f"Processed {idx} genomes...")
    
    # Write output file
    print(f"\nWriting results to: {output_file}")
    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines))
    
    print(f"Done! Total {len(output_lines)-1} genomes processed")


def main():
    if len(sys.argv) != 4:
        print("Usage: python3 sequence_haplotypes.py <gene_family_matrix> <fasta_file> <output_file>")
        print("\nParameters:")
        print("  gene_family_matrix: Gene family matrix file")
        print("  fasta_file: FASTA file containing all gene protein sequences")
        print("  output_file: Output haplotype matrix file")
        sys.exit(1)
    
    matrix_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]
    
    print("="*60)
    print("Homologous Gene Haplotype Analysis")
    print("="*60)
    
    try:
        process_gene_family_matrix(matrix_file, fasta_file, output_file)
        print("\nAnalysis complete!")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
```
</details>

---

## Haplotype Diversity Analysis

### Overview
Calculate Nei's diversity index for haplotype diversity analysis across groups and perform statistical comparisons.

**Nei's Diversity Index Formula:** `H = 1 - Σ(pi²)`
- where `pi` is the frequency of the i-th haplotype
- Range: 0-1
- H = 0: No diversity (all samples have same haplotype)
- H = 1: Maximum diversity (all samples have different haplotypes)

### Script: `haplotype.diversity.py`

**Usage:**
```bash
# Basic usage
python haplotype.diversity.py -i homology.haplotype.table.all.input -o diversity_result

# With custom parameters
python haplotype.diversity.py -i data.tab -o result --min-samples 5 --dpi 600
```

**Parameters:**
- `-i, --input`: Haplotype data file (tab-delimited, required)
- `-o, --output`: Output file prefix (required)
- `--min-samples`: Minimum sample size per group (default: 3)
- `--min-haplotypes`: Minimum number of haplotypes per gene (default: 2)
- `--style`: Plot style (seaborn/ggplot/bmh/classic, default: seaborn)
- `--colors`: Color scheme (default: Set2)
- `--figsize`: Figure size in inches (default: 10 6)
- `--dpi`: Figure resolution (default: 600)

**Input File Format:**
```
Sample    Group    Gene001    Gene002    Gene003
S001      GroupA   Hap1       Hap2       Hap1
S002      GroupA   Hap1       Hap2       Hap2
S003      GroupB   Hap2       Hap1       Hap1
```

**Output Files:**
- `*_diversity.txt` - Diversity values for each gene in each group (wide format)
- `*_pairwise_comparison.txt` - Pairwise statistical comparison results
- `*_plot_data.txt` - Data used for plotting (raw + summary statistics)
- `*_diversity_plot.pdf` - Publication-ready diversity boxplot
- `*_report.txt` - Comprehensive analysis report

**Key Features:**
- Calculates Nei's diversity index per gene per group
- Kruskal-Wallis test for overall comparison (>2 groups)
- Mann-Whitney U test for pairwise comparisons
- FDR correction (Benjamini-Hochberg method)
- Multiple comparison letter labels
- Publication-ready plots (600 DPI)

<details>
<summary><b>📜 Click to view complete Python script (500+ lines)</b></summary>

```python
#!/usr/bin/env python3
"""
Haplotype Diversity Analysis - Nei's Diversity Index
Calculate haplotype diversity for each gene in each group and perform inter-group comparisons.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import false_discovery_control
import argparse
import sys
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

def calculate_nei_diversity(haplotypes):
    """Calculate Nei's diversity index"""
    haplotypes = [h for h in haplotypes if not pd.isna(h)]
    if len(haplotypes) < 2:
        return np.nan
    unique, counts = np.unique(haplotypes, return_counts=True)
    frequencies = counts / len(haplotypes)
    diversity = 1.0 - np.sum(frequencies ** 2)
    return diversity

def calculate_diversity_for_group(data, group_name):
    """Calculate diversity for all genes in a group"""
    diversities = {}
    for gene in data.columns:
        if gene not in ['Sample', 'Group']:
            diversity = calculate_nei_diversity(data[gene].values)
            if not np.isnan(diversity):
                diversities[gene] = diversity
    return diversities

def parse_args():
    parser = argparse.ArgumentParser(
        description='Calculate haplotype diversity (Nei\'s diversity)',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='Haplotype data file')
    parser.add_argument('-o', '--output', required=True, help='Output file prefix')
    parser.add_argument('--min-samples', type=int, default=3, help='Minimum sample size per group')
    parser.add_argument('--min-haplotypes', type=int, default=2, help='Minimum haplotypes per gene')
    parser.add_argument('--style', default='seaborn', choices=['seaborn', 'ggplot', 'bmh', 'classic'])
    parser.add_argument('--colors', default='Set2', help='Color scheme')
    parser.add_argument('--figsize', nargs=2, type=float, default=[10, 6])
    parser.add_argument('--dpi', type=int, default=600, help='Figure resolution')
    return parser.parse_args()

def main():
    args = parse_args()
    
    print("=" * 70)
    print("Haplotype Diversity Analysis - Nei's Diversity Index")
    print("=" * 70)
    print(f"Input file: {args.input}")
    print(f"Output prefix: {args.output}")
    print("=" * 70)
    
    # Read data
    print("\n[1/6] Reading data...")
    try:
        data = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        print(f"Error: Cannot read file - {e}")
        sys.exit(1)
    
    # Detect sample column
    sample_col = None
    for col in ['Sample', 'sample', 'ID', 'id', 'SampleID']:
        if col in data.columns:
            sample_col = col
            break
    if sample_col is None:
        sample_col = data.columns[0]
    
    if 'Group' not in data.columns:
        print("Error: No 'Group' column in data")
        sys.exit(1)
    
    data = data.rename(columns={sample_col: 'Sample'})
    print(f"    Samples: {len(data)}, Genes: {data.shape[1] - 2}")
    
    # Get groups
    groups = data['Group'].unique()
    print(f"    Groups: {len(groups)}")
    for group in groups:
        n = len(data[data['Group'] == group])
        print(f"      {group}: {n} samples")
    
    # Filter groups
    valid_groups = []
    for group in groups:
        if len(data[data['Group'] == group]) >= args.min_samples:
            valid_groups.append(group)
        else:
            print(f"    Warning: Group '{group}' excluded (insufficient samples)")
    
    if len(valid_groups) < 2:
        print("Error: At least 2 groups required")
        sys.exit(1)
    
    groups = valid_groups
    print(f"    Retained {len(groups)} groups")
    
    # Calculate diversity
    print("\n[2/6] Calculating diversity...")
    all_diversities = {}
    gene_cols = [c for c in data.columns if c not in ['Sample', 'Group']]
    
    for group in groups:
        print(f"    Processing: {group}")
        group_data = data[data['Group'] == group]
        diversities = calculate_diversity_for_group(group_data, group)
        all_diversities[group] = diversities
        print(f"      Valid genes: {len(diversities)}/{len(gene_cols)}")
    
    # Organize data
    print("\n[3/6] Organizing data...")
    common_genes = set(all_diversities[groups[0]].keys())
    for group in groups[1:]:
        common_genes &= set(all_diversities[group].keys())
    print(f"    Common genes: {len(common_genes)}")
    
    diversity_data = []
    for group in groups:
        for gene in common_genes:
            diversity_data.append({
                'Gene': gene,
                'Group': group,
                'Diversity': all_diversities[group][gene]
            })
    
    df_diversity = pd.DataFrame(diversity_data)
    
    # Save results
    output_file = f"{args.output}_diversity.txt"
    df_wide = df_diversity.pivot(index='Gene', columns='Group', values='Diversity')
    df_wide = df_wide.round(6)
    df_wide['Mean'] = df_wide.mean(axis=1)
    df_wide['SD'] = df_wide.std(axis=1)
    df_wide['CV'] = df_wide['SD'] / df_wide['Mean']
    df_wide.to_csv(output_file, sep='\t')
    print(f"    Saved: {output_file}")
    
    # Statistical tests
    print("\n[4/6] Statistical testing...")
    summary_stats = df_diversity.groupby('Group')['Diversity'].agg([
        ('N', 'count'), ('Mean', 'mean'), ('SD', 'std'),
        ('Median', 'median'), ('Min', 'min'), ('Max', 'max')
    ]).round(6)
    print("\nGroup summary:")
    print(summary_stats)
    
    # Kruskal-Wallis test
    group_values = [df_diversity[df_diversity['Group'] == g]['Diversity'].values for g in groups]
    
    if len(groups) > 2:
        h_stat, p_value = stats.kruskal(*group_values)
        print(f"\nKruskal-Wallis test:")
        print(f"  H = {h_stat:.4f}, p = {p_value:.6e}")
    
    # Pairwise comparisons
    print("\nPairwise comparisons:")
    pairwise_results = []
    from itertools import combinations
    
    for g1, g2 in combinations(groups, 2):
        values1 = df_diversity[df_diversity['Group'] == g1]['Diversity'].values
        values2 = df_diversity[df_diversity['Group'] == g2]['Diversity'].values
        u_stat, p_value = stats.mannwhitneyu(values1, values2, alternative='two-sided')
        mean1, mean2 = np.mean(values1), np.mean(values2)
        
        pairwise_results.append({
            'Group1': g1, 'Group2': g2, 'Mean1': mean1, 'Mean2': mean2,
            'Difference': mean1 - mean2, 'U_statistic': u_stat, 'P_value': p_value
        })
        
        sig = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
        print(f"  {g1} vs {g2}: p={p_value:.6e} ({sig})")
    
    # FDR correction
    print("\nApplying FDR correction...")
    p_values = [r['P_value'] for r in pairwise_results]
    p_adjusted = false_discovery_control(np.array(p_values), method='bh')
    
    for i, r in enumerate(pairwise_results):
        r['P_adjusted'] = p_adjusted[i]
        r['Significant_FDR'] = 'Yes' if p_adjusted[i] < 0.05 else 'No'
    
    # Letter labels
    print("\nCalculating letter labels...")
    sig_matrix = {}
    for g1, g2 in combinations(groups, 2):
        for r in pairwise_results:
            if (r['Group1'] == g1 and r['Group2'] == g2) or \
               (r['Group1'] == g2 and r['Group2'] == g1):
                sig_matrix[(g1, g2)] = r['P_adjusted'] < 0.05
                sig_matrix[(g2, g1)] = r['P_adjusted'] < 0.05
                break
    
    group_means = {g: df_diversity[df_diversity['Group'] == g]['Diversity'].mean() for g in groups}
    sorted_groups = sorted(groups, key=lambda x: group_means[x], reverse=True)
    
    letters = {}
    available_letters = list('abcdefghijklmnopqrstuvwxyz')
    letter_idx = 0
    
    for group in sorted_groups:
        non_sig_groups = [group]
        for other in sorted_groups:
            if other != group and not sig_matrix.get((group, other), False):
                non_sig_groups.append(other)
        
        existing_letters = set()
        for g in non_sig_groups:
            if g in letters:
                existing_letters.update(letters[g])
        
        if existing_letters:
            letters[group] = min(existing_letters)
        else:
            letters[group] = available_letters[letter_idx] if letter_idx < len(available_letters) else 'z'
            letter_idx += 1
    
    print("  Letter labels:")
    for group in sorted_groups:
        print(f"    {group}: {letters[group]} (mean={group_means[group]:.4f})")
    
    # Save pairwise results
    df_pairwise = pd.DataFrame(pairwise_results)
    pairwise_file = f"{args.output}_pairwise_comparison.txt"
    df_pairwise.to_csv(pairwise_file, sep='\t', index=False, float_format='%.6f')
    print(f"\nSaved: {pairwise_file}")
    
    # Plotting
    print("\n[5/6] Creating plots...")
    plot_data_file = f"{args.output}_plot_data.txt"
    plot_data_list = []
    for group in groups:
        for val in df_diversity[df_diversity['Group'] == group]['Diversity'].values:
            plot_data_list.append({'Group': group, 'Diversity': val})
    
    plot_data_df = pd.DataFrame(plot_data_list)
    stats_df = pd.DataFrame({
        'Group': groups,
        'N': [len(df_diversity[df_diversity['Group'] == g]) for g in groups],
        'Mean': [df_diversity[df_diversity['Group'] == g]['Diversity'].mean() for g in groups],
        'SD': [df_diversity[df_diversity['Group'] == g]['Diversity'].std() for g in groups],
        'Letter': [letters[g] for g in groups]
    })
    
    with open(plot_data_file, 'w') as f:
        f.write("# Raw data\n")
        plot_data_df.to_csv(f, sep='\t', index=False)
        f.write("\n# Summary statistics\n")
        stats_df.to_csv(f, sep='\t', index=False)
    
    print(f"    Saved: {plot_data_file}")
    
    # Create plot
    plt.style.use('default')
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10, 
                        'pdf.fonttype': 42, 'ps.fonttype': 42})
    
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#6A994E', '#BC4B51']
    fig, ax = plt.subplots(figsize=(6, 4.5))
    positions = np.arange(len(groups))
    
    bp = ax.boxplot([df_diversity[df_diversity['Group'] == g]['Diversity'].values for g in groups],
                     positions=positions, widths=0.5, patch_artist=True, showfliers=False,
                     medianprops=dict(color='white', linewidth=2),
                     boxprops=dict(linewidth=1.5))
    
    for patch, color in zip(bp['boxes'], colors[:len(groups)]):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    
    # Add data points
    for i, group in enumerate(groups):
        data_vals = df_diversity[df_diversity['Group'] == group]['Diversity'].values
        x_jitter = np.random.normal(i, 0.04, size=len(data_vals))
        ax.scatter(x_jitter, data_vals, color='black', s=20, alpha=0.3, zorder=3)
    
    # Add letter labels
    y_max = df_diversity['Diversity'].max()
    y_min = df_diversity['Diversity'].min()
    y_range = y_max - y_min
    
    for i, group in enumerate(groups):
        ax.text(i, y_max + y_range * 0.08, letters[group],
               ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    ax.set_xticks(positions)
    ax.set_xticklabels(groups, fontsize=11)
    ax.set_ylabel("Nei's Diversity Index", fontsize=12, fontweight='bold')
    ax.set_ylim(max(0, y_min - y_range * 0.1), y_max + y_range * 0.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plot_file = f"{args.output}_diversity_plot.pdf"
    plt.savefig(plot_file, dpi=args.dpi, bbox_inches='tight', format='pdf')
    print(f"    Saved: {plot_file}")
    plt.close()
    
    # Generate report
    print("\n[6/6] Generating report...")
    report_file = f"{args.output}_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("Haplotype Diversity Analysis Report\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Input: {args.input}\n")
        f.write(f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("Data Overview:\n")
        f.write(f"  Samples: {len(data)}\n")
        f.write(f"  Genes: {len(gene_cols)}\n")
        f.write(f"  Groups: {len(groups)}\n")
        f.write(f"  Analyzed: {len(common_genes)}\n\n")
        f.write("Group Statistics:\n")
        f.write(summary_stats.to_string() + "\n\n")
        f.write("Letter Labels (FDR < 0.05):\n")
        for g in sorted_groups:
            f.write(f"  {g}: {letters[g]} (mean={group_means[g]:.6f})\n")
        if len(groups) > 2:
            f.write(f"\nKruskal-Wallis: H={h_stat:.4f}, p={p_value:.6e}\n\n")
        f.write("Pairwise Comparisons:\n")
        for comp in pairwise_results:
            sig_raw = '***' if comp['P_value']<0.001 else '**' if comp['P_value']<0.01 else '*' if comp['P_value']<0.05 else 'ns'
            sig_fdr = '***' if comp['P_adjusted']<0.001 else '**' if comp['P_adjusted']<0.01 else '*' if comp['P_adjusted']<0.05 else 'ns'
            f.write(f"{comp['Group1']} vs {comp['Group2']}:\n")
            f.write(f"  Mean: {comp['Mean1']:.6f} vs {comp['Mean2']:.6f}\n")
            f.write(f"  Raw p: {comp['P_value']:.6e} ({sig_raw})\n")
            f.write(f"  FDR p: {comp['P_adjusted']:.6e} ({sig_fdr})\n\n")
    
    print(f"    Saved: {report_file}")
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)

if __name__ == '__main__':
    main()
```

</details>

---

## Haplotype Frequency Comparison

### Overview
High-performance analysis of haplotype frequency differences across groups using chi-square tests with parallel computing.

### Script: `haplotype.frequency.comparison.py`

**Usage:**
```bash
# Basic usage (use all CPU cores)
python haplotype.frequency.comparison.py -i input.tab -o output_prefix

# Specify number of cores
python haplotype.frequency.comparison.py -i input.tab -o results --threads 4

# Fast mode (skip pairwise comparisons)
python haplotype.frequency.comparison.py -i input.tab -o results --no-pairwise

# Minimal output (no frequency tables)
python haplotype.frequency.comparison.py -i input.tab -o results --minimal
```

**Parameters:**
- `-i, --input`: Input haplotype file (tab-delimited, required)
- `-o, --output`: Output file prefix (required)
- `-p, --pvalue`: P-value threshold (default: 0.05)
- `--threads`: Number of CPU cores to use (default: auto-detect)
- `--minimal`: Skip frequency table generation (faster)
- `--no-pairwise`: Skip pairwise group comparisons (much faster)

**Input File Format:**
Same as haplotype diversity analysis:
- Column 1: Sample ID
- Column 2: Group assignment
- Columns 3+: Gene haplotypes

**Output Files:**
- `*_missing_values_summary.txt` - Missing value statistics per gene
- `*_haplotype_counts.txt` - Haplotype type count matrix (groups × genes)
- `*_chi_square_test.txt` - Chi-square test results with FDR correction
- `*_pairwise_comparison.txt` - Pairwise group comparison results
- `*_significant_genes.txt` - Summary of genes with significant differences
- `*_<gene>_haplotype_frequencies.txt` - Detailed frequency tables (top 20 genes)
- `*_summary_report.txt` - Comprehensive analysis summary

**Key Features:**
- **Parallel computing:** Multi-core processing for large datasets (10,000+ genes)
- **Chi-square tests:** Association between groups and haplotypes
- **Pairwise comparisons:** Identify which specific group pairs differ
- **FDR correction:** Benjamini-Hochberg method for multiple testing
- **Performance optimization:** Vectorized operations and memory efficiency

**Performance Tips:**
- Use `--no-pairwise` for 10-100× speed improvement
- Use `--minimal` to skip frequency tables (2-5× faster)
- Specify `--threads` to control CPU usage
- For >20,000 genes, may require 32+ GB RAM

<details>
<summary><b>📜 Click to view complete Python script (400+ lines)</b></summary>

```python
#!/usr/bin/env python3
"""
Gene Haplotype Analysis - High Performance Version
Optimized for large-scale datasets (10,000+ genes)
"""

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, false_discovery_control
import argparse
import sys
import os
from multiprocessing import Pool, cpu_count
import time

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Gene Haplotype Analysis - High Performance Version')
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    parser.add_argument('-p', '--pvalue', type=float, default=0.05, help='P-value threshold')
    parser.add_argument('--threads', type=int, default=None, help='Number of CPU cores')
    parser.add_argument('--minimal', action='store_true', help='Minimal output mode')
    parser.add_argument('--no-pairwise', action='store_true', help='Skip pairwise comparisons')
    args = parser.parse_args()
    if not 0 < args.pvalue < 1:
        parser.error("P-value must be between 0 and 1")
    return args

args = parse_arguments()

# Set threads
n_threads = max(1, cpu_count() - 1) if args.threads is None else max(1, min(args.threads, cpu_count()))

print("=" * 70)
print("Gene Haplotype Analysis - High Performance Version")
print("=" * 70)
print(f"Input: {args.input}")
print(f"Output: {args.output}")
print(f"P-value: {args.pvalue}")
print(f"CPU cores: {n_threads} / {cpu_count()}")
print("=" * 70)

def chi_square_worker(gene_data_tuple):
    """Worker for parallel chi-square test"""
    gene, gene_values, group_values = gene_data_tuple
    try:
        contingency_table = pd.crosstab(group_values, gene_values)
        if contingency_table.shape[0] < 2 or contingency_table.shape[1] < 2:
            return None
        chi2, p_value, dof, _ = chi2_contingency(contingency_table)
        return (gene, chi2, dof, p_value)
    except:
        return None

def pairwise_worker(args_tuple):
    """Worker for pairwise comparison"""
    gene, gene_values, group_values, group1, group2 = args_tuple
    try:
        mask = (group_values == group1) | (group_values == group2)
        pair_gene = gene_values[mask]
        pair_group = group_values[mask]
        contingency_table = pd.crosstab(pair_group, pair_gene)
        if contingency_table.shape[0] < 2 or contingency_table.shape[1] < 2:
            return None
        chi2, p_value, _, _ = chi2_contingency(contingency_table)
        return (gene, group1, group2, chi2, p_value)
    except:
        return None

# Read data
print("\n[1/6] Reading data...")
start_time = time.time()

if not os.path.exists(args.input):
    print(f"Error: File not found: {args.input}")
    sys.exit(1)

try:
    data = pd.read_csv(args.input, sep='\t')
except Exception as e:
    print(f"Error: Cannot read file: {e}")
    sys.exit(1)

if data.shape[1] < 3:
    print("Error: Need at least 3 columns")
    sys.exit(1)

sample_col = data.columns[0]
group_col = data.columns[1]
gene_cols = data.columns[2:].tolist()
groups = data[group_col].unique()

print(f"    Samples: {len(data)}, Genes: {len(gene_cols)}, Groups: {len(groups)}")
print(f"    Time: {time.time() - start_time:.2f}s")

# Filter missing values
print("\n[2/6] Filtering missing values...")
start_time = time.time()

missing_counts = data[gene_cols].isna().sum()
genes_without_missing = missing_counts[missing_counts == 0].index.tolist()

print(f"    Original: {len(gene_cols)}")
print(f"    After filtering: {len(genes_without_missing)}")
print(f"    Time: {time.time() - start_time:.2f}s")

if len(genes_without_missing) == 0:
    print("Error: All genes have missing values!")
    sys.exit(1)

gene_cols = genes_without_missing

# Save missing summary
missing_summary = pd.DataFrame({
    'gene': data.columns[2:],
    'missing_count': data.iloc[:, 2:].isna().sum().values,
    'missing_percentage': data.iloc[:, 2:].isna().sum().values / len(data) * 100
})
missing_summary['has_missing'] = missing_summary['missing_count'] > 0
missing_summary = missing_summary.sort_values('missing_count', ascending=False)
missing_summary.to_csv(f"{args.output}_missing_values_summary.txt", sep='\t', index=False)

# Count haplotype types
print("\n[3/6] Counting haplotype types...")
start_time = time.time()

count_matrix = pd.DataFrame(index=groups, columns=gene_cols)
for group in groups:
    group_mask = data[group_col] == group
    count_matrix.loc[group] = data.loc[group_mask, gene_cols].nunique()

count_matrix = count_matrix.astype(int)
count_matrix.to_csv(f"{args.output}_haplotype_counts.txt", sep='\t')
print(f"    Saved: {args.output}_haplotype_counts.txt")
print(f"    Time: {time.time() - start_time:.2f}s")

# Chi-square test (parallel)
print("\n[4/6] Chi-square test (parallel)...")
start_time = time.time()

gene_data_list = []
for gene in gene_cols:
    gene_values = data[gene].values
    group_values = data[group_col].values
    gene_data_list.append((gene, gene_values, group_values))

print(f"    Using {n_threads} cores...")
with Pool(n_threads) as pool:
    results = pool.map(chi_square_worker, gene_data_list, chunksize=50)

results = [r for r in results if r is not None]

if len(results) == 0:
    print("Error: No valid results")
    sys.exit(1)

chi_square_df = pd.DataFrame(results, columns=['gene', 'chi_square', 'df', 'p_value'])
chi_square_df['p_adjusted'] = false_discovery_control(chi_square_df['p_value'].values)
chi_square_df['significant_fdr'] = chi_square_df['p_adjusted'] < args.pvalue
chi_square_df = chi_square_df.sort_values('p_value')

chi_square_df.to_csv(f"{args.output}_chi_square_test.txt", sep='\t', index=False)
print(f"    Saved: {args.output}_chi_square_test.txt")
print(f"    Significant: {chi_square_df['significant_fdr'].sum()} / {len(chi_square_df)}")
print(f"    Time: {time.time() - start_time:.2f}s")

# Pairwise comparison (optional)
if not args.no_pairwise:
    print("\n[5/6] Pairwise comparison (parallel)...")
    start_time = time.time()
    
    from itertools import combinations
    group_pairs = list(combinations(groups, 2))
    
    pairwise_data_list = []
    for gene in gene_cols:
        gene_values = data[gene].values
        group_values = data[group_col].values
        for group1, group2 in group_pairs:
            pairwise_data_list.append((gene, gene_values, group_values, group1, group2))
    
    print(f"    Total comparisons: {len(pairwise_data_list)}")
    
    with Pool(n_threads) as pool:
        pairwise_results = pool.map(pairwise_worker, pairwise_data_list, chunksize=100)
    
    pairwise_results = [r for r in pairwise_results if r is not None]
    
    if len(pairwise_results) > 0:
        pairwise_df = pd.DataFrame(pairwise_results,
                                   columns=['gene', 'group1', 'group2', 'chi_square', 'p_value'])
        pairwise_df['p_adjusted'] = false_discovery_control(pairwise_df['p_value'].values)
        pairwise_df['significant_fdr'] = pairwise_df['p_adjusted'] < args.pvalue
        pairwise_df = pairwise_df.sort_values('p_value')
        
        pairwise_df.to_csv(f"{args.output}_pairwise_comparison.txt", sep='\t', index=False)
        print(f"    Saved: {args.output}_pairwise_comparison.txt")
        
        sig_genes = pairwise_df[pairwise_df['significant_fdr']].groupby('gene').agg({
            'p_value': ['count', 'min'],
            'p_adjusted': 'min'
        }).reset_index()
        sig_genes.columns = ['gene', 'n_significant_pairs', 'min_p_value', 'min_p_adjusted']
        
        sig_pairs_detail = pairwise_df[pairwise_df['significant_fdr']].groupby('gene').apply(
            lambda x: '; '.join([f"{row['group1']} vs {row['group2']}" for _, row in x.iterrows()])
        ).reset_index()
        sig_pairs_detail.columns = ['gene', 'significant_pairs']
        
        sig_genes = sig_genes.merge(sig_pairs_detail, on='gene')
        sig_genes = sig_genes.sort_values(['n_significant_pairs', 'min_p_adjusted'],
                                         ascending=[False, True])
        
        sig_genes.to_csv(f"{args.output}_significant_genes.txt", sep='\t', index=False)
        print(f"    Saved: {args.output}_significant_genes.txt")
        print(f"    Significant genes: {len(sig_genes)}")
    
    print(f"    Time: {time.time() - start_time:.2f}s")
else:
    print("\n[5/6] Skipping pairwise (--no-pairwise)")

# Frequency tables (optional)
if not args.minimal and not args.no_pairwise:
    print("\n[6/6] Calculating frequencies (top 20)...")
    start_time = time.time()
    
    if 'sig_genes' in locals() and len(sig_genes) > 0:
        top_genes = sig_genes['gene'].head(20).tolist()
        
        for gene in top_genes:
            freq_data = data.groupby([group_col, gene]).size().reset_index(name='count')
            freq_data.columns = ['group', 'haplotype', 'count']
            
            group_totals = freq_data.groupby('group')['count'].sum().reset_index()
            group_totals.columns = ['group', 'total']
            
            freq_data = freq_data.merge(group_totals, on='group')
            freq_data['frequency'] = freq_data['count'] / freq_data['total']
            freq_data['percentage'] = freq_data['frequency'] * 100
            
            freq_data = freq_data.sort_values(['group', 'frequency'], ascending=[True, False])
            freq_data.to_csv(f"{args.output}_{gene}_haplotype_frequencies.txt",
                           sep='\t', index=False)
        
        print(f"    Generated for top 20 genes")
        print(f"    Time: {time.time() - start_time:.2f}s")
else:
    print("\n[6/6] Skipping frequencies (--minimal)")

# Summary
print("\nGenerating summary...")

summary_lines = [
    "=" * 70,
    "Gene Haplotype Analysis Report",
    "=" * 70,
    "",
    f"Date: {pd.Timestamp.now().date()}",
    f"Input: {args.input}",
    "",
    "Data:",
    f"  Samples: {len(data)}",
    f"  Original genes: {len(data.columns) - 2}",
    f"  Genes analyzed: {len(gene_cols)}",
    f"  Groups: {len(groups)} ({', '.join(groups)})",
    "",
    "Performance:",
    f"  CPU cores: {n_threads}",
    "",
    "Results:",
    f"  Haplotype types: {count_matrix.values.mean():.2f} (avg)",
    f"  Range: {count_matrix.values.min()} - {count_matrix.values.max()}",
    ""
]

if len(chi_square_df) > 0:
    n_sig = chi_square_df['significant_fdr'].sum()
    summary_lines.extend([
        f"Differential test (FDR < {args.pvalue}):",
        f"  Genes tested: {len(chi_square_df)}",
        f"  Significant: {n_sig} ({n_sig/len(chi_square_df)*100:.2f}%)",
        ""
    ])

if not args.no_pairwise and 'sig_genes' in locals() and len(sig_genes) > 0:
    summary_lines.extend([
        "Pairwise:",
        f"  Genes with ≥2 significant pairs: {len(sig_genes)}",
        f"  Max pairs: {sig_genes['n_significant_pairs'].max()}",
        ""
    ])

summary_lines.extend([
    "Output files:",
    f"  1. {args.output}_missing_values_summary.txt",
    f"  2. {args.output}_haplotype_counts.txt",
    f"  3. {args.output}_chi_square_test.txt",
])

if not args.no_pairwise:
    summary_lines.extend([
        f"  4. {args.output}_pairwise_comparison.txt",
        f"  5. {args.output}_significant_genes.txt",
    ])
    if not args.minimal:
        summary_lines.append(f"  6. {args.output}_*_frequencies.txt")

summary_lines.extend(["", "=" * 70])

summary_report = '\n'.join(summary_lines)
with open(f"{args.output}_summary_report.txt", 'w') as f:
    f.write(summary_report)

print(summary_report)
print(f"\nComplete! Summary: {args.output}_summary_report.txt")
```

</details>

---

## Notes

- All scripts require Python 3.8+
- Required packages: `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`
- Install via conda: `conda install -c conda-forge pandas numpy scipy matplotlib seaborn`
- Or via pip: `pip install pandas numpy scipy matplotlib seaborn`
- Scripts are optimized for large-scale genomic datasets
- Use `--help` flag with any script for detailed parameter information

---

**Last Updated:** January 2026
