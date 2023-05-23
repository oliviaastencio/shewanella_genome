#! /usr/bin/env bash

## COMPARE GENOMES
. ~soft_bio_267/initializes/init_autoflow  

project_path=`pwd`
data_path=$project_path'/data'
script_path=$project_path'/script'
genome_analysis_path=$SCRATCH/genome_analysis

GENOME_FILES=$data_path/genomes_problem   #"/mnt/home/users/pab_001_uma/oliastencio/micro_lab3/Pdp11_genome_analysis/TarSynFlow/pathogenic_strains/genomes"
PROTEOME=$genome_analysis_path/genes_identification/Tarsynflow/proteome/prots_clean.fasta   #"/mnt/home/users/pab_001_uma/oliastencio/micro_lab3/Pdp11_genome_analysis/TarSynFlow/pathogenic_strains/proteome/prots_clean.fasta"
while read QUERY
do
	sed 's/>/>Q_/g' "$GENOME_FILES/$QUERY" > "$GENOME_FILES/e_$QUERY"
	while read REF
	do
		OUTPUT=$genome_analysis_path/genes_identification/Tarsynflow/comps/$QUERY/$REF
		mkdir -p $OUTPUT
		sed 's/>/>R_/g' "$GENOME_FILES/$REF" > "$GENOME_FILES/e_$REF"
		VARS=`echo "genomes_path=$GENOME_FILES
		query_genome=e_$QUERY
		ref_genome=e_$REF
		proteins_dataset=$PROTEOME
		coverage_identity_query=85+85
		coverage_identity_ref=85.0+85.0
		genomic=''" | sed 's/^/$/' | tr "\n" "," | tr -d "\t"`
		echo $VARS
		AutoFlow -w $script_path/workflow_genome_sinteny_with_proteins -c 10 -s -n 'cal' -t '1-00:00:00' $1 -V $VARS -o $OUTPUT
	done < $data_path/gen_refs
done < $data_path/gen_queries
