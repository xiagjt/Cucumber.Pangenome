```markdown
# Genomic Repeat Analysis Pipeline

This document provides the commands used for building a database and running repeat analysis on the genomic data.

## Software Installation

Before running the commands, ensure that you have the required software installed. You can use Conda to install RepeatModeler and RepeatMasker as follows:

```bash
# Create a new conda environment
conda create -n repeat_analysis python=3.8

# Activate the environment
conda activate repeat_analysis

# Install RepeatModeler
conda install -c bioconda repeatmodeler

# Install RepeatMasker
conda install -c bioconda repeatmasker
```

## Commands

### 1. Build the Database

To create a database from the genomic sequence file, use the following command:

```bash
BuildDatabase -name CG104 CG104.Nuclear.Add.fa
```

### 2. Run RepeatModeler

Next, run RepeatModeler to identify and classify repeats in the genomic data. The command below specifies the necessary parameters and directories:

```bash
RepeatModeler -pa 60 -LTRStruct -ninja_dir /home/GuanJianTao/Software/NINJA-0.95-cluster_only/NINJA -rmblast_dir /home/GuanJianTao/micromamba/envs/repeatmodeler/bin/ -trf_dir /home/GuanJianTao/micromamba/envs/repeatmodeler/bin/ -database CG104
```

### 3. Run RepeatMasker

Finally, apply RepeatMasker to mask the identified repeats using the library you created. The command is as follows:

```bash
RepeatMasker -lib CG104-families.fa -pa 30 -html -gff -dir Masker.out CG104.Nuclear.Add.fa
```

## Requirements

- **Conda**: Ensure that Conda is installed on your system.
- **RepeatModeler**: Installed via Conda, with paths to necessary directories specified.
- **RepeatMasker**: Installed via Conda, with access to the repeat library (`CG104-families.fa`).

## Output

The output files will be generated in the specified directories, with HTML and GFF formats provided for easy viewing and further analysis.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

```

You can further adjust any sections as necessary to fit your project's needs!
