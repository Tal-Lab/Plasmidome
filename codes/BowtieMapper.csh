#!/bin/tcsh
#$ -o ../logs/Output_Bowtie2
#$ -e ../logs/Error_Bowtie2

cd $3
bowtie2 -k 2 -x $4 -1 $1 -2 $2 -S $5
echo "$6_double.bam"
samtools view -bS $5 > "$6_double.bam"
samtools sort "$6_double.bam" -o "$6_double.sorted.bam"
bcftools mpileup -Ou $7 "$6_double.sorted.bam" | bcftools call -Ou -mv - > "$6_double.raw.bcf"
bcftools view "$6_double.raw.bcf"
cd ../../codes