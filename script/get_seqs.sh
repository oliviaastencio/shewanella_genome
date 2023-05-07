#! /usr/bin/env bash

mkdir -p genes_identification
mkdir -p genes_identification/Tarsynflow
mkdir -p genes_identification/Tarsynflow/proteome
#mkdir genomes
wget  'https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=Shewanella' -O  genes_identification/Tarsynflow/proteome/prots.fasta
#gunzip proteome/prots.fasta.gz
