#!/bin/tcsh
#$ -o $6/Output_BWA
#$ -e $6/Error_BWA

recycler_path=$RECYCLER

python $recycler_path/make_fasta_from_fastg.py -g $1 -o $2
bwa/bwa index $2
bwa/bwa mem $2 $4 $5 | samtools view -buS - > $3/reads_pe.bam
samtools view -bF 0x0800 $3/reads_pe.bam > $3/reads_pe_primary.bam
samtools sort $3/reads_pe_primary.bam -o $3/reads_pe_primary.sort.bam
samtools index $3/reads_pe_primary.sort.bam
