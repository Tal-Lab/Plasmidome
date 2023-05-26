#!/bin/tcsh
#$ -o ../logs/plsdb_blast_output
#$ -e ../logs/plsdb_blast_error

blast_path=$BLAST
plsdb_link=$PLSDB

mkdir -p DBs/PLSDB/reference
cd DBs/PLSDB/reference

wget $plsdb_link -O plsdb.fna

$blast_path/makeblastdb -in plsdb.fna -dbtype 'nucl' -out plsdb -logfile 'PLSDB_log.log'

cd ../../../Output

$blast_path/blastn -query "/filtered_plasmids.fasta" -db "../DBs/PLSDB/reference/plsdb" -evalue 0.001 -out "plsdb.csv" -outfmt '6 qseqid sseqid stitle evalue length pident mismatch score qcov qstart qend sstart send qseq sseq'

