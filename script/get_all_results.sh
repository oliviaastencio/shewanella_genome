#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --time=02:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

## COMPARE GENOMES
#source ~soft_cvi_114/initializes/init_autoflow
. ~soft_bio_267/initializes/init_autoflow
###################################################################
## FUNCTIONS
##################################################################
function get_annotations { #$1 => protein_list $2 => min_count
	prot_list=`awk '{
	if ($2 >= '$2')
		print $0;
	}' 'genes_identification/Tarsynflow/results/'$1'_freq_proteins' | cut -f 1 | sed 's/$/+OR+/g' | tr -d "\n" |sed 's/+OR+$//g'`
	wget "https://rest.uniprot.org/uniprotkb/stream?=true&fields=accession,length,gene_names,reviewed,protein_name,ec,organism_name,go&format=tsv&query=accession:$prot_list" -O 'genes_identification/Tarsynflow/results/'$1'_specific_filtered_annotated_prots'
} 

function get_filtered_list { #$1 => protein_list $2 => min_count $3 => output list
        awk '{
        if ($2 >= '$2')
                print $0;
        }' 'genes_identification/Tarsynflow/results/'$1'_freq_proteins' | cut -f 1 > 'genes_identification/Tarsynflow/results/'$3'_filtered_prots'
}

function get_frequency { #$1 => input_file $2 => protein_list $3 => freq_output_file 
	while read prot
	do
		rm temp_stats
		while read REF
		do
			grep -c $prot $comparation_folder/$QUERY/$REF'/circos_0000/'$1 >> temp_stats
		done < data/gen_refs_1
		echo -e "$prot\t`grep -v -c -w '0' temp_stats`" >> 'genes_identification/Tarsynflow/results/'$3
	done < 'genes_identification/Tarsynflow/results/protein_lists/'$2
	rm temp_stats
}

function get_protein_list { #$1 => input_file, $2 temp_out_file $3 => final_protein_list
        rm $comparation_folder/$QUERY/$2
        while read REF
        do
                cut -f 1 -d ' ' $comparation_folder/$QUERY/$REF/circos_0000/$1 >> $comparation_folder/$QUERY/$2
        done < data/gen_refs_1
        sort -u  $comparation_folder/$QUERY/$2 > genes_identification/Tarsynflow/results/$3
}

#################################################################
## MAIN
################################################################### 
comparation_folder=genes_identification/Tarsynflow/comps
min_count=5

rm -rf genes_identification/Tarsynflow/results
mkdir genes_identification/Tarsynflow/results
mkdir genes_identification/Tarsynflow/results/images
mkdir genes_identification/Tarsynflow/results/matches_analysis
while read QUERY
do
        while read REF
        do	
		# CHECK LAUNCHED WORKFLOWS
		ex_path=$comparation_folder/$QUERY/$REF
		echo $ex_path
		flow_logger -e $ex_path -r all #-l -w
		# COPY CIRCOS IMAGES
		circos_path="$comparation_folder/$QUERY/$REF/circos_0000"
		compared_genomes=`cat $circos_path/compared_genomes`
		cp $circos_path/circos.png genes_identification/Tarsynflow/results/images/$compared_genomes".png"
		cut -f 1 -d ' ' $circos_path"/specific_seqs_REF.txt" | sed "s/^/$REF /g" >> genes_identification/Tarsynflow/results/matches_analysis/pairs_genome_prot
		cut -f 1 -d ' ' $circos_path"/specific_seqs_QUERY.txt" | sed "s/^/$QUERY /g" >> genes_identification/Tarsynflow/results/matches_analysis/pairs_genome_prot
        done < data/gen_refs_1

	# GET PROTEIN LISTS
	mkdir genes_identification/Tarsynflow/results/protein_lists
	get_protein_list specific_seqs_QUERY.txt query_specific_matches protein_lists/Query_common_specific_proteins
	get_protein_list specific_seqs_REF.txt ref_specific_matches protein_lists/Ref_common_specific_proteins
	get_protein_list shared_seqs_SURE.txt sure_shared_matches protein_lists/sure_common_shared_proteins
	get_protein_list shared_seqs_PUTATIVE.txt putative_shared_matches protein_lists/putative_common_shared_proteins
	get_protein_list not_matching_SURE.txt sure_not_matching_matches protein_lists/sure_common_not_matching_proteins
	get_protein_list not_matching_PUTATIVE.txt putative_not_matching_matches protein_lists/putative_common_not_matching_proteins

	# GET PROTEIN FREQUENCY LIST
	get_frequency specific_seqs_QUERY.txt Query_common_specific_proteins Query_freq_proteins
	get_frequency specific_seqs_REF.txt Ref_common_specific_proteins Ref_freq_proteins
	get_frequency shared_seqs_SURE.txt sure_common_shared_proteins sure_shared_freq_proteins
	get_frequency shared_seqs_PUTATIVE.txt putative_common_shared_proteins putative_shared_freq_proteins
	get_frequency not_matching_SURE.txt sure_common_not_matching_proteins sure_not_matching_freq_proteins
	get_frequency not_matching_PUTATIVE.txt putative_common_not_matching_proteins putative_not_matching_freq_proteins

	#FILTER PROTEIN LIST BY FREQ AND RETRIEVE ANNOTATIONS	
	get_annotations Query $min_count
	get_annotations Ref $min_count
	get_filtered_list sure_shared $min_count sure_shared
	get_filtered_list putative_shared $min_count putative_shared
	get_filtered_list sure_not_matching $min_count sure_not_matching
	get_filtered_list putative_not_matching $min_count putative_not_matching
	
done < data/gen_queries_1

. ~soft_bio_267/initializes/init_ruby
script/scripts_Tarsynflow/pairs2matrix.rb genes_identification/Tarsynflow/results/matches_analysis/pairs_genome_prot genes_identification/Tarsynflow/results/matches_analysis/matrix_genome_prot
module load R
script/scripts_Tarsynflow/plot_heatmap.R -f genes_identification/Tarsynflow/results/matches_analysis/matrix_genome_prot -o genes_identification/Tarsynflow/results/matches_analysis/heatmap.pdf