# Orthogroup Classification and Statistical Analysis

## Overview
Classify orthogroups based on presence-absence patterns across different sample groups and sources (Alt vs Ref), and perform comprehensive statistical analysis.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Input Files](#input-files)
- [Analysis Workflow](#analysis-workflow)
- [Output Files](#output-files)
- [Example Results](#example-results)
- [Use Cases](#use-cases)

---

## Installation

### Requirements
- Python 3.8+
- pandas

### Install Dependencies
```bash
pip install pandas
```

---

## Usage

### Basic Command
```bash
python orthogroup_analysis.py \
  -a Syntenic.groups.GeneCount.tsv \
  -b Sample.group.xls \
  -o Out
```

### Parameters

| Parameter | Short | Required | Description |
|-----------|-------|----------|-------------|
| `--count` | `-a` | Yes | Gene count matrix file |
| `--sample` | `-b` | Yes | Sample classification information file |
| `--output` | `-o` | Yes | Output file prefix |

### Help
```bash
python orthogroup_analysis.py --help
```

---

## Input Files

### 1. Gene Count Matrix (`-a`)

Tab-delimited file with orthogroups as rows and samples as columns.

**Format:**
```
Orthogroup    Sample1    Sample2    Sample3    Sample4
SG1     2          1          0          3
SG2     0          0          1          1
SG3     5          4          6          5
SG4     1          1          1          1
```

**Structure:**
- **First column**: Orthogroup IDs (e.g., OG0000001, OG0000002)
- **Other columns**: Gene copy numbers per sample
- **Values**: 
  - ≥ 1 indicates presence of genes in that sample
  - 0 indicates absence

**Notes:**
- This is typically the output from OrthoFinder (`Orthogroups.GeneCount.tsv`)
- Sample names in header must match those in the sample classification file
- Missing samples will generate warnings but won't stop analysis

---

### 2. Sample Classification File (`-b`)

Tab-delimited file with three required columns defining sample metadata.

**Format:**
```
Sample     Group      Source
Sample1    Wild       Alt
Sample2    Wild       Alt
Sample3    Cultivated Ref
Sample4    Cultivated Ref
Sample5    EastAsian  Alt
Sample6    Eurasian   Ref
```

**Required Columns:**

#### Column 1: Sample
- Sample identifier (must match column names in gene count matrix)
- Examples: "Sample1", "Accession_001", "Genome_ABC"
- Must be unique
- No spaces allowed (will be automatically stripped)

#### Column 2: Group
- Ecological or geographical group classification
- Examples: "Wild", "Cultivated", "EastAsian", "Eurasian", "GroupA", "GroupB"
- Used for within-source group analysis
- Can have multiple groups per source

#### Column 3: Source
- Data source classification
- **Only two values allowed:**
  - `Alt`: Alternative/New genome assemblies (your newly sequenced genomes)
  - `Ref`: Reference genome assemblies (published reference genomes)
- This distinction is used for Alt vs Ref comparison

**Example with Real Data:**
```
Sample              Group        Source
CG104_Alt          EastAsian    Alt
CG37_Alt           Eurasian     Alt
9930_Ref           EastAsian    Ref
Gy14_Ref           Eurasian     Ref
Wild_CN01          Wild         Alt
Cultivar_USA       Cultivated   Ref
```

**Notes:**
- File can be Excel (.xls, .xlsx) or tab-delimited text (.tsv, .txt)
- Column names are flexible but must include these three categories
- Spaces in values will be automatically stripped
- Case-sensitive for Group and Source values

---

## Analysis Workflow

### Step 1: Alt vs Ref Classification

Classify orthogroups based on presence in Alt vs Ref samples.

**Categories:**

| Category | Definition | Example Use |
|----------|------------|-------------|
| **Alt_specific** | Present in ≥1 Alt sample, absent in all Ref samples | Genes unique to new assemblies |
| **Ref_specific** | Present in ≥1 Ref sample, absent in all Alt samples | Genes only in reference genomes |
| **Shared** | Present in both Alt and Ref samples | Core genes common to both |
| **Absent_both** | Absent in both Alt and Ref samples | Should be rare (quality check) |

**Logic:**
- Presence defined as gene count ≥ 1
- Absence defined as gene count = 0
- "Present in Alt" means present in at least one Alt sample
- "Present in Ref" means present in at least one Ref sample

---

### Step 2: Ecological Group Analysis (Alt Samples)

Analyze orthogroup distribution across ecological groups within Alt samples only.

**Categories:**

| Category | Definition | Example |
|----------|------------|---------|
| **All_groups_shared** | Present in all groups within Alt | Core genes across all Alt groups |
| **Group_specific** | Present only in one specific group | Wild_specific, Cultivated_specific |
| **Group_shared** | Present in specific combinations of groups | Wild_Cultivated_shared |
| **All_groups_absent** | Absent in all groups | Rare, quality check |

**Example:**
If Alt samples have three groups (Wild, Cultivated, EastAsian):
- **Wild_specific**: Present in Wild samples only
- **Cultivated_specific**: Present in Cultivated samples only
- **Wild_Cultivated_shared**: Present in Wild and Cultivated, absent in EastAsian
- **All_groups_shared**: Present in all three groups

---

### Step 3: Ecological Group Analysis (Ref Samples)

Same analysis as Step 2, but for Ref samples only.

This allows comparison of gene distribution patterns between Alt and Ref datasets.

---

### Step 4: Ecological Group Analysis (All Samples)

Same analysis combining both Alt and Ref samples together.

This provides an overall view of gene distribution across all ecological groups regardless of source.

---

### Step 5: Cross-tabulation Analysis

Generate contingency tables showing relationships between:
- Alt/Ref categories (from Step 1)
- Ecological group categories (from Steps 2-4)

This reveals patterns such as:
- How many Alt-specific genes are Wild-specific?
- How many shared genes are present across all groups?
- Distribution of Ref-specific genes across ecological groups

---

## Output Files

### Summary Files

#### 1. Alt vs Ref Summary
**File:** `{output}_Alt_vs_Ref_summary.tsv`

```
Category          Count
Alt_specific      1234
Ref_specific      567
Shared            8901
Absent_both       45
```

#### 2. Alt Group Summary
**File:** `{output}_Alt_Group_summary.tsv`

```
Category                       Count
All_groups_shared             5678
All_groups_absent             123
Wild_specific                 234
Cultivated_specific           345
Wild_Cultivated_shared        456
EastAsian_Cultivated_shared   89
```

#### 3. Ref Group Summary
**File:** `{output}_Ref_Group_summary.tsv`

Same format as Alt Group Summary, but for Ref samples.

#### 4. All Group Summary
**File:** `{output}_All_Group_summary.tsv`

Same format, but combining Alt and Ref samples.

---

### Orthogroup Lists

#### Alt vs Ref Lists
- `{output}_Alt_specific_OGs.txt` - One orthogroup ID per line
- `{output}_Ref_specific_OGs.txt` - One orthogroup ID per line
- `{output}_Alt_Ref_shared_OGs.txt` - One orthogroup ID per line

**Example:**
```
OG0000001
OG0000015
OG0000234
OG0001567
```

#### Group-Specific Lists (Alt)
- `{output}_Alt_All_groups_shared_OGs.txt`
- `{output}_Alt_{Group}_specific_OGs.txt` (e.g., `Out_Alt_Wild_specific_OGs.txt`)
- `{output}_Alt_{Groups}_shared_OGs.txt` (e.g., `Out_Alt_Wild_Cultivated_shared_OGs.txt`)

#### Group-Specific Lists (Ref)
- `{output}_Ref_All_groups_shared_OGs.txt`
- `{output}_Ref_{Group}_specific_OGs.txt`
- `{output}_Ref_{Groups}_shared_OGs.txt`

#### Group-Specific Lists (All)
- `{output}_All_All_groups_shared_OGs.txt`
- `{output}_All_{Group}_specific_OGs.txt`
- `{output}_All_{Groups}_shared_OGs.txt`

---

### Detailed Classification Table

**File:** `{output}_OG_classification.tsv`

Complete classification of every orthogroup across all three analysis levels.

```
Orthogroup    Alt_Ref_Category    Alt_Group_Category       Ref_Group_Category       All_Group_Category
OG0000001     Alt_specific        Wild_specific            All_groups_absent        Wild_specific
OG0000002     Shared              All_groups_shared        All_groups_shared        All_groups_shared
OG0000003     Ref_specific        All_groups_absent        Cultivated_specific      Cultivated_specific
OG0000004     Shared              Wild_specific            Eurasian_specific        Wild_Eurasian_shared
```

**Columns:**
- **Orthogroup**: Orthogroup ID
- **Alt_Ref_Category**: Classification from Alt vs Ref analysis
- **Alt_Group_Category**: Classification from Alt group analysis
- **Ref_Group_Category**: Classification from Ref group analysis
- **All_Group_Category**: Classification from All group analysis

**Use Cases:**
- Filter orthogroups by specific criteria
- Identify patterns across classification schemes
- Input for downstream enrichment analysis
- Visualization of gene distribution

---

### Cross-tabulation Tables

#### 1. Alt Groups Cross-tab
**File:** `{output}_cross_tab_Alt_groups.tsv`

```
Alt_Ref_Category    Wild_specific    Cultivated_specific    All_groups_shared    ...
Alt_specific        234              156                    0                    ...
Ref_specific        0                0                      0                    ...
Shared              12               45                     5678                 ...
```

Shows distribution of Alt/Ref categories across Alt ecological groups.

#### 2. Ref Groups Cross-tab
**File:** `{output}_cross_tab_Ref_groups.tsv`

Same format, for Ref ecological groups.

#### 3. All Groups Cross-tab
**File:** `{output}_cross_tab_All_groups.tsv`

Same format, for All ecological groups combined.

**Use Cases:**
- Identify patterns between source and ecological classification
- Quality control (e.g., Alt_specific should not appear in Ref groups)
- Quantify overlap between different classification schemes

---

## Example Results

### Complete Analysis Example

**Input:**
- 15,000 orthogroups
- 50 Alt samples (20 Wild, 15 Cultivated, 15 EastAsian)
- 30 Ref samples (10 Wild, 10 Cultivated, 10 Eurasian)

**Output Summary:**

```
Alt vs Ref Classification:
  Alt_specific:    1,234 (8.2%)
  Ref_specific:    567 (3.8%)
  Shared:          13,156 (87.7%)
  Absent_both:     43 (0.3%)

Alt Group Classification:
  All_groups_shared:  8,901 (59.3%)
  Wild_specific:      234 (1.6%)
  Cultivated_specific: 345 (2.3%)
  EastAsian_specific:  189 (1.3%)
  Wild_Cultivated_shared: 456 (3.0%)
  [Additional combinations...]

Ref Group Classification:
  All_groups_shared:  9,234 (61.6%)
  Wild_specific:      123 (0.8%)
  Cultivated_specific: 198 (1.3%)
  Eurasian_specific:  287 (1.9%)
  [Additional combinations...]
```

---

## Use Cases

### 1. Pan-genome Analysis
**Goal:** Identify core, accessory, and private genes

**Approach:**
- **Core genes**: Use "All_groups_shared" from All analysis
- **Accessory genes**: Use group-shared categories
- **Private genes**: Use group-specific categories

**Example:**
```bash
# Core genes present in all samples
cat Out_All_All_groups_shared_OGs.txt

# Genes specific to wild accessions
cat Out_All_Wild_specific_OGs.txt
```

---

### 2. Comparative Genomics
**Goal:** Compare gene content between Alt (new) and Ref (published) genomes

**Approach:**
- Identify Alt-specific genes (potentially novel genes)
- Identify Ref-specific genes (potentially missed in Alt)
- Quantify shared gene content

**Example:**
```bash
# Novel genes in new assemblies
cat Out_Alt_specific_OGs.txt

# Genes potentially missed in new assemblies
cat Out_Ref_specific_OGs.txt
```

---

### 3. Population Structure Analysis
**Goal:** Analyze gene distribution across ecological/geographical groups

**Approach:**
- Use group-specific and group-shared categories
- Compare patterns between Alt and Ref datasets
- Use cross-tabulation to identify associations

**Example:**
```bash
# Genes unique to East Asian group
cat Out_Alt_EastAsian_specific_OGs.txt

# Genes shared between Wild and Cultivated only
cat Out_Alt_Wild_Cultivated_shared_OGs.txt
```

---

### 4. Domestication Studies
**Goal:** Identify genes lost or gained during domestication

**Approach:**
- Compare Wild vs Cultivated groups
- Use Alt analysis to focus on newly sequenced genomes
- Identify Wild-specific and Cultivated-specific genes

**Example:**
```bash
# Genes lost during domestication
cat Out_Alt_Wild_specific_OGs.txt

# Genes gained during domestication
cat Out_Alt_Cultivated_specific_OGs.txt
```

---

### 5. Quality Control
**Goal:** Assess genome assembly completeness

**Approach:**
- Check "Absent_both" category (should be minimal)
- Compare Alt vs Ref shared genes (should be high)
- Identify systematic differences between sources

**Quality Indicators:**
- **Good**: Shared > 85%, Absent_both < 1%
- **Acceptable**: Shared > 75%, Absent_both < 2%
- **Needs review**: Shared < 75% or Absent_both > 2%

---

## Key Features

### Presence-based Classification
- Gene count ≥ 1 = present
- Gene count = 0 = absent
- Robust to copy number variation

### Multi-level Analysis
- Alt vs Ref comparison
- Within-Alt group analysis
- Within-Ref group analysis
- Combined All analysis
- Independent but complementary views

### Comprehensive Output
- Summary statistics for quick overview
- Detailed orthogroup lists for downstream analysis
- Classification table for filtering and visualization
- Cross-tabulation for pattern discovery

### Group Combinations
- Automatically identifies all possible group combinations
- Useful for identifying genes with specific distribution patterns
- Examples: genes in Wild+Cultivated but not EastAsian

### Quality Checks
- Warns about sample mismatches between input files
- Reports number of samples per group
- Identifies potentially problematic categories (e.g., high Absent_both)

---

## Troubleshooting

### Common Issues

#### 1. Sample Name Mismatch
**Error:** Warning about missing samples

**Solution:**
- Ensure sample names in count matrix header match exactly with Sample column in classification file
- Check for extra spaces, tabs, or special characters
- Use `head -1 count.tsv` and compare with `cut -f1 sample.txt`

#### 2. No Alt or Ref Samples
**Error:** Analysis fails if all samples are one source

**Solution:**
- Ensure you have both Alt and Ref samples in classification file
- Check Source column values (must be exactly "Alt" or "Ref")

#### 3. Only One Group
**Warning:** Limited group analysis

**Solution:**
- Need at least 2 groups for meaningful group analysis
- Group-specific analysis will work but group-shared will be empty

#### 4. High Absent_both Count
**Concern:** Many orthogroups absent in all samples

**Possible Causes:**
- Low-quality gene predictions
- Stringent filtering in OrthoFinder
- Orthogroups from outgroup species

**Solution:**
- Review OrthoFinder parameters
- Check if outgroup species were included
- Consider filtering these orthogroups from analysis

---

## Notes

### Important Considerations

- **Gene counts**: Only presence/absence matters, not copy number
- **Group definitions**: User-defined, flexible for different research questions
- **Source labels**: Must be exactly "Alt" and "Ref" (case-sensitive)
- **Sample matching**: Critical for correct analysis
- **Output volume**: Can be large with many groups and combinations

### Best Practices

1. **Validate inputs**: Check sample name consistency before running
2. **Document groups**: Keep record of what each group represents
3. **Start simple**: Begin with major groups, then refine
4. **Use cross-tabs**: Helpful for understanding relationships
5. **Filter results**: Focus on categories relevant to your question

### Performance

- **Speed**: Fast for typical datasets (<30 seconds for 15K orthogroups)
- **Memory**: Low memory requirements (< 1GB for 15K orthogroups)
- **Scalability**: Can handle 50K+ orthogroups without issues

---

## Citation

If you use this script in your research, please cite:

```
Zhang et al. (2025) Pangenome-resolved structural variation drives adaptation 
and trait evolution in cucumber. Nature Genetics (in press). 
Manuscript ID: NG-A69677R
```

---

## Contact

For questions, issues, or suggestions:
- Open an issue on GitHub: https://github.com/xiagjt/Cucumber.Pangenome/issues
- Contact the corresponding author (see publication)

---

**Last Updated:** January 2025
**Version:** 1.0
