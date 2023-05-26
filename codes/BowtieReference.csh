#!/bin/tcsh
#$ -o ../logs/Output_Bowtie2
#$ -e ../logs/Error_Bowtie2

bowtie2-build -f $1 $2


