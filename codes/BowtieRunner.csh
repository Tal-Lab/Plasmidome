#!/bin/tcsh
#$ -o ../logs/Output_Bowtie2
#$ -e ../logs/Error_Bowtie2

cp "../res/filtered_plasmids.fasta" "../Output/plasmids_copy.fasta"
echo 'Creating fasta with concatenated plasmid candidates'
seqkit concat "../res/filtered_plasmids.fasta" "../Output/plasmids_copy.fasta" > "../Output/plasmids_double.fasta"
mkdir "../Output/reference"
cp "../Output/plasmids_double.fasta" "../Output/reference/plasmids_double.fasta"
python "/bowtie_runner.py" "../Output/plasmids_double.fasta"


