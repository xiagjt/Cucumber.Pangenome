```markdown
# NLR-Annotator and BRAKER2 Pipeline

This document outlines the commands used for NLR annotation and gene prediction using BRAKER2.

## NLR-Annotator

To annotate NLR genes, use the following command. This command reads IDs from `IDs.xls` and processes each line:

```bash
cat IDs.xls | while read line
do
    java -jar /home/GuanJianTao/Software/NLR-Annotator/NLR-Annotator-v2.1b.jar -i ${line}.final.fasta -x /home/GuanJianTao/Software/NLR-Annotator/src/mot.txt -y /home/GuanJianTao/Software/NLR-Annotator/src/store.txt -o ${line}.mask.genome.fa.output.txt -t 30
done
```

## BRAKER2

For gene prediction using BRAKER2, set the environment variables and run the following command:

```bash
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
braker.pl --cores=60 --genome=CG104.20K.fa --prot_seq=journal.pbio.3001124.s013.fasta --epmode --gff3 --min_contig=10000 --workingdir=CG104.braker.out
```

## Requirements

- **Java**: Ensure that Java is installed and properly configured to run NLR-Annotator.
- **BRAKER2**: Make sure that BRAKER2 and its dependencies are installed.
- **Input Files**: Ensure that the input files (`IDs.xls`, `${line}.final.fasta`, `CG104.20K.fa`, and `journal.pbio.3001124.s013.fasta`) are available in the working directory.

## Output

The output files will be generated as specified in the commands, including annotated NLR outputs and gene predictions in GFF3 format.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

```
