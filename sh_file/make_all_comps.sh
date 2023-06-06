#! /usr/bin/env bash


## COMPARE GENOMES
. ~soft_bio_267/initializes/init_autoflow  

data_path=$1     #$project_path'/data'
script_path=$2    #$project_path'/script'
genome_analysis_path=$3 ##$SCRATCH/genome_analysis

GENOME_FILES=$data_path/genomes_problem   
PROTEOME=$genome_analysis_path/genes_identification/Tarsynflow/proteome/prots_clean.fasta   
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
		AutoFlow -w $script_path/workflow_genome_sinteny_with_proteins -c 10 -s -n 'cal' -t '1-00:00:00' $1 -V $VARS -o $OUTPUT
	done < $data_path/gen_refs
done < $data_path/gen_queries
