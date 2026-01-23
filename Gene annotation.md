# Genome Annotation Pipeline

This workflow describes the structural annotation of a genome using a combination of transcriptome-based evidence, *de novo* prediction, and protein homology.

## 1. Environment Setup

We use `Conda` to manage the various software suites required for this pipeline.

```bash
# Create a master environment for annotation tools
conda create -n annotation_env -c bioconda -c conda-forge \
    hisat2 stringtie trinity pasa braker2 augustus \
    genomethreader evidencemodeler samtools gffread bedtools bioawk -y

# Activate the environment
conda activate annotation_env

```

---

## 2. Transcriptome Evidence Generation

### 2.1 Reference-Based Assembly

This approach aligns RNA-seq reads to the genome to identify transcribed regions.

```bash
# 1. Build Index
hisat2-build Csa406_genome_np2_r2.fa.soft_masked.fa Csa406_index

# 2. Alignment
hisat2 -x Csa406_index -1 $fq1 -2 $fq2 -p 10 --dta 2>&1 | tee $stat -S $out

# 3. Sort and Convert to BAM
samtools view -bS $out | samtools sort -@ 20 -o $bam -

# 4. Assemble Transcripts
stringtie $bam -p 20 -o $out_gtf

# 5. Merge assemblies (TACO or StringTie --merge)
taco_run gtf_files.txt

```

### 2.2 Reference-Free (De novo) Assembly

Useful for capturing transcripts that may be difficult to align due to genome assembly gaps.

```bash
# Assemble using Trinity
Trinity --seqType fq --left $left --right $right --CPU 50 \
        --max_memory 50G --normalize_reads --full_cleanup \
        --output trinity_out

```

### 2.3 Refinement with PASA

PASA aligns the assembled transcripts to the genome to model gene structures.

```bash
# Run PASA pipeline (for both ref-based and de novo outputs)
Launch_PASA_pipeline.pl -c alignAssembly.config -C -R --ALIGNER gmap \
    -g Csa406_genome_np2_r2.fa.soft_masked.fa -t $input_fasta --CPU 40

```

---

## 3. Protein and Ab Initio Prediction

### 3.1 BRAKER2 (Augustus)

BRAKER2 uses RNA-seq (BAM) and protein evidence to train Augustus automatically.

```bash
export OPENBLAS_NUM_THREADS=1

braker.pl --genome=$genome --prot_seq=$protein_db --bam=$bam \
          --softmasking --cores 60 --etpmode

# Convert output to GFF
gffread braker/augustus.hints.gtf -o- > augustus.hints.gff

```

### 3.2 Homology-based (GenomeThreader)

Aligns known protein sequences to the target genome.

```bash
# Alignment
gth -genomic $genome -protein $protein_parts -intermediate -xmlout \
    -prseedlength 20 -gcmincoverage 80 -prminmatchlen 20 -o gth_output.gz

# Generate Consensus
gthconsensus -gff3out -intermediate -o Broad.consensus.gff *.gz

```

---

## 4. Evidence Integration (EvidenceModeler)

EvidenceModeler (EVM) combines all previous outputs into a single, weighted consensus gene set.

### 4.1 Preparation

Convert all various formats (Braker, PASA, GenomeThreader) into EVM-compatible GFF3.

```bash
# Convert Braker/Augustus
perl braker_GTF_to_EVM_GFF3.pl braker.gtf > braker.EVM.gff
perl augustus_GFF3_to_EVM_GFF3.pl augustus.hints.gff > augustus.EVM.gff

# Convert GenomeThreader
perl genomeThreader_to_evm_gff3.pl Broad.consensus.homo.final.gff > Broad.concensus.final.EVM.gff

```

### 4.2 Running EVM

Partition the genome, run the model, and recombine.

```bash
# 1. Partition Genome
partition_EVM_inputs.pl --genome $genome --gene_predictions Prediction.gff \
    --transcript_alignments $rna_gff --segmentSize 1000000 --overlapSize 100000

# 2. Write and Execute Commands
write_EVM_commands.pl --genome $genome --weights weights.txt \
    --gene_predictions Prediction.gff --transcript_alignments $rna_gff \
    --output_file_name evm.out --partitions partitions_list.out > commands.list

parallel --j 80 < commands.list

# 3. Recombine and Convert to GFF3
recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out --genome $genome

```

---

## 5. Post-Processing and Filtering

Final step to remove low-quality annotations (e.g., genes encoding proteins shorter than 50 amino acids).

```bash
# Extract CDS sequences
gffread EVM.all.gff -g $genome -y tr_cds.fa

# Identify short proteins (< 50 aa)
bioawk -c fastx '$seq < 50 {print $comment}' tr_cds.fa | cut -d '=' -f 2 > short_list.txt

# Filter the final GFF
grep -v -w -f short_list.txt EVM.all.gff > final_annotation.filter.gff

```
