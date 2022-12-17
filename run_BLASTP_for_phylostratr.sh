#!/bin/bash

focal_fasta=$1
uniprot_seqs=$2
threads=$3

### MAIN ###

for uniprot_fasta in $(find $uniprot_seqs -type f); do
    tag=$(basename $uniprot_fasta)
    echo "Comparison with $uniprot_fasta"
    mkdir ${tag}_blast_out
    cd ${tag}_blast_out
    nohup makeblastdb -dbtype prot -in $uniprot_fasta -out $tag
    wait
    echo -e 'qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tscore' > ${tag}.tab
    nohup blastp -query $focal_fasta -db $tag -num_threads $threads -outfmt '6 qseqid sseqid qstart qend sstart send evalue score' >> ${tag}.tab
    wait
    cd ..
done