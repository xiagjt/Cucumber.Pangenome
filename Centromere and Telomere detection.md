This pipeline is designed for the high-throughput identification of **Centromeric Regions** across a pangenome scale (130 samples) using sequence alignment, distance-based merging, and length-based filtering.

---

# Centromere Identification Pipeline

This workflow automates the process of building BLAST databases for multiple genome assemblies, performing sequence alignment with centromere-specific repeats, and identifying the "core" centromere region on each chromosome.

## 1. Environment Setup

The pipeline requires **NCBI-BLAST+** for alignment and **Bedtools** for genomic interval manipulation.

```bash
# Create and activate the environment
conda create -n centromere_env -c bioconda blast bedtools -y
conda activate centromere_env

```

---

## 2. BLAST Database Construction

We iterate through all 130 genome assemblies in the pangenome directory to create individual nucleotide BLAST databases.

<details>
<summary>Click to expand: make_blastdb.sh</summary>

```bash
#!/usr/bin/env bash

# Directories
GENOME_DIR="/home/Public.Resources/Genomes/PanGenome.S130.Final.dir"
DB_DIR="/vol2/ZhangYiFan2023/workspace/centromere/blastdatabase"

for sample_dir in "$GENOME_DIR"/*; do
    sample=$(basename "$sample_dir")
    fasta_file="$sample_dir/genome/${sample}.final.new.fasta"
    out_dir="$DB_DIR/$sample"

    if [[ ! -f "$fasta_file" ]]; then
        echo "⚠️ Skipping $sample: $fasta_file not found"
        continue
    fi

    mkdir -p "$out_dir"
    
    # Generate BLAST database
    makeblastdb -in "$fasta_file" \
                -dbtype nucl \
                -out "$out_dir/$sample" \
                -parse_seqids \
                -title "$sample"
done

```

</details>

---

## 3. High-Throughput BLAST Alignment

Align the centromere query sequence against all generated databases using `blastn`.

<details>
<summary>Click to expand: run_blast.sh</summary>

```bash
#!/usr/bin/env bash

# Settings
QUERY="/vol2/ZhangYiFan2023/workspace/centromere/centromere.fa"
DB_DIR="/vol2/ZhangYiFan2023/workspace/centromere/blastdatabase"
OUT_DIR="/vol2/ZhangYiFan2023/workspace/centromere/blast_result"
mkdir -p "$OUT_DIR"

for sample_dir in "$DB_DIR"/*; do
    sample=$(basename "$sample_dir")
    db_prefix="$sample_dir/$sample"

    if [[ ! -f "${db_prefix}.nhr" ]]; then
        echo "Skipping $sample: Database not found"
        continue
    fi
    
    out_file="$OUT_DIR/${sample}.blast.txt"

    # Execute blastn with tabular output format (outfmt 6)
    blastn \
        -query "$QUERY" \
        -db "$db_prefix" \
        -out "$out_file" \
        -outfmt 6 \
        -num_threads 40
done

```

</details>

---

## 4. Result Filtering and BED Conversion

Filter the BLAST hits based on E-value, alignment length, and specific repeat types. Results are converted to BED format for downstream proximity analysis.

<details>
<summary>Click to expand: blast_filter.sh</summary>

```bash
#!/usr/bin/env bash

BLAST_DIR="/vol2/ZhangYiFan2023/workspace/centromere/blast_result"
OUT_DIR="/vol2/ZhangYiFan2023/workspace/centromere/bam.dir"
mkdir -p "$OUT_DIR"

for file in "$BLAST_DIR"/*.txt; do
    sample=$(basename "$file" .blast.txt)
    outfile="$OUT_DIR/${sample}.bed"

    # Filter: type_III repeats, E-value < 1e-5, and alignment length > 100bp
    # Adjust coordinates to ensure Start < End for BED format
    awk '{if ($1=="type_III" && $11 < 1e-5 && $4 > 100) {
            OFS="\t";
            if ($9 < $10) {
                print $2, $9, $10
            } else {
                print $2, $10, $9
            }
        }
    }' "$file" | sort -k1,1 -k2,2n > "$outfile"
done

```

</details>

---

## 5. Hit Merging and Core Region Selection

Centromeres are often composed of clusters of repeats. We merge hits within **100 kb** of each other and identify the longest merged block as the candidate centromere for each chromosome.

<details>
<summary>Click to expand: blast_merge.sh</summary>

```bash
#!/usr/bin/env bash

indir="/vol2/ZhangYiFan2023/workspace/centromere/bam.dir"
merged_dir="/vol2/ZhangYiFan2023/workspace/centromere/bam.dir/merged"
mkdir -p "$merged_dir"

output_all="All_Germplasms_Merged_MaxRegion.tsv"
echo -e "Germplasm\tChr\tStart\tEnd\tLength" > "$output_all"

for bed in "$indir"/*.bed; do
    germ=$(basename "$bed" .bed)
    merged_bed="$merged_dir/${germ}.merged.bed"

    # Merge intervals within 100,000 bp distance
    bedtools merge -i "$bed" -d 100000 > "$merged_bed"

    # Select the longest merged region per chromosome (putative centromere core)
    awk '
        {
            len=$3-$2
            key=$1
            if(len > maxlen[key]) {
                maxlen[key]=len
                line[key]=$0
            }
        }
        END {
            for(k in line) print line[k]
        }
    ' "$merged_bed" | sort -k1,1 > "$merged_dir/${germ}.max_region.bed"

    # Compile final summary
    awk -v g="$germ" '{
        len=$3-$2;
        print g"\t"$1"\t"$2"\t"$3"\t"len
    }' "$merged_dir/${germ}.max_region.bed" >> "$output_all"
done

```

</details>

---

# Telomere Identification Pipeline

https://github.com/jamiemcg/TelomereSearch

