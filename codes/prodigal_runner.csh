#!/bin/tcsh
#$ -o ../logs/Output_Prodigal
#$ -e ../logs/Error_Prodigal

prodigal -i "../Output/plasmids_double.fasta" -o "../Output/plasmids_double_proteins.gbk" -a "../Output/plasmids_double_proteins.faa" -p anon