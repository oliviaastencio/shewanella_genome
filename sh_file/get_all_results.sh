#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
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
	}' $1'_freq_proteins' | cut -f 1 | sed 's/$/+OR+/g' | tr -d "\n" |sed 's/+OR+$//g'`
	wget "https://rest.uniprot.org/uniprotkb/stream?=true&fields=accession,length,gene_names,reviewed,protein_name,ec,organism_name,go&format=tsv&query=accession:$prot_list" -O $1'_specific_filtered_annotated_prots'
} 

function get_filtered_list { #$1 => protein_list $2 => min_count $3 => output list
        awk '{
        if ($2 >= '$2')
                print $0;
        }' $1'_freq_proteins' | cut -f 1 > $3'_filtered_prots'
}

function get_frequency { #$1 => input_file $2 => protein_list $3 => freq_output_file 
	while read prot
	do
		rm temp_stats
		while read REF
		do
			grep -c $prot $comparation_folder/$QUERY/$REF'/circos_0000/'$1 >> temp_stats
		done < $data/gen_refs
		echo -e "$prot\t`grep -v -c -w '0' temp_stats`" >> $3
	done < $2
	rm temp_stats
}

function get_protein_list { #$1 => input_file, $2 temp_out_file $3 => final_protein_list
        rm $comparation_folder/$QUERY/$2
        while read REF
        do
                cut -f 1 -d ' ' $comparation_folder/$QUERY/$REF/circos_0000/$1 >> $comparation_folder/$QUERY/$2
        done < $data/gen_refs
        sort -u  $comparation_folder/$QUERY/$2 > $3
}

#################################################################
## MAIN
################################################################### 
output=$1
data=$2 
comparation_folder=$3  #genes_identification/Tarsynflow/comps
script=$4

min_count=5

if [ "$5" == "get_protein" ]; then  

rm -rf $output
mkdir $output
mkdir $output/images
mkdir $output/matches_analysis

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
		cp $circos_path/circos.png $output/images/$compared_genomes".png"
		cut -f 1 -d ' ' $circos_path"/specific_seqs_REF.txt" | sed "s/^/$REF /g" >> $output/matches_analysis/pairs_genome_prot
		cut -f 1 -d ' ' $circos_path"/specific_seqs_QUERY.txt" | sed "s/^/$QUERY /g" >> $output/matches_analysis/pairs_genome_prot
		
    done < $data/gen_refs

    GET PROTEIN LISTS
	mkdir -p $output/protein_lists
	get_protein_list specific_seqs_QUERY.txt query_specific_matches  $output/protein_lists/Query_common_specific_proteins
	get_protein_list specific_seqs_REF.txt ref_specific_matches  $output/protein_lists/Ref_common_specific_proteins
	get_protein_list shared_seqs_SURE.txt sure_shared_matches  $output/protein_lists/sure_common_shared_proteins
	get_protein_list shared_seqs_PUTATIVE.txt putative_shared_matches  $output/protein_lists/putative_common_shared_proteins
	get_protein_list not_matching_SURE.txt sure_not_matching_matches  $output/protein_lists/sure_common_not_matching_proteins
	get_protein_list not_matching_PUTATIVE.txt putative_not_matching_matches  $output/protein_lists/putative_common_not_matching_proteins

	# GET PROTEIN FREQUENCY LIST
	get_frequency specific_seqs_QUERY.txt $output/protein_lists/Query_common_specific_proteins $output/Query_freq_proteins
	get_frequency specific_seqs_REF.txt $output/protein_lists/Ref_common_specific_proteins $output/Ref_freq_proteins
	get_frequency shared_seqs_SURE.txt $output/protein_lists/sure_common_shared_proteins $output/sure_shared_freq_proteins
	get_frequency shared_seqs_PUTATIVE.txt $output/protein_lists/putative_common_shared_proteins $output/putative_shared_freq_proteins
	get_frequency not_matching_SURE.txt $output/protein_lists/sure_common_not_matching_proteins $output/sure_not_matching_freq_proteins
	get_frequency not_matching_PUTATIVE.txt $output/protein_lists/putative_common_not_matching_proteins $output/putative_not_matching_freq_proteins

done < $data/gen_queries
fi


if [ "$5" == "get_annotations" ]; then 
while read QUERY
do

	#FILTER PROTEIN LIST BY FREQ AND RETRIEVE ANNOTATIONS	
	get_annotations $output/Query $min_count
	get_annotations $output/Ref $min_count
	get_filtered_list $output/sure_shared $min_count $output/sure_shared
	get_filtered_list $output/putative_shared $min_count $output/putative_shared
	get_filtered_list $output/sure_not_matching $min_count $output/sure_not_matching
	get_filtered_list $output/putative_not_matching $min_count $output/putative_not_matching
	
done < $data/gen_queries

. ~soft_bio_267/initializes/init_ruby
$script/scripts_Tarsynflow/pairs2matrix.rb $output/matches_analysis/pairs_genome_prot $output/matches_analysis/matrix_genome_prot
module load R
$script/scripts_Tarsynflow/plot_heatmap.R -f $output/matches_analysis/matrix_genome_prot -o $output/matches_analysis/heatmap.pdf
fi 

###Nota:si la descarga de $1'_specific_filtered_annotated_prots da error, ejecutar sin enviar al sistema de colas. 