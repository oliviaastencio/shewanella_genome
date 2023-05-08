#! /usr/bin/env bash

######################################
#####FUNCTIONS
#####################################

##################################################
######Shewanella complete genomes database 
#################################################
function genomes_database {
	output=$1
	module load ruby/2.4.1 
	rm -r $output/genomes_Shewanella
	mkdir -p $output/genomes_Shewanella/sequence_download
	
	main_folder=`pwd` 

	prokaryotes_genomes_csv=$output/prokaryotes.csv  #all shewanella genomes from https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/shewanella
	input=$output/genomes_Shewanella/sequence_download #sequence_download
	rm $output/genomes_link
	cat $prokaryotes_genomes_csv | tr ',' '\t'| cut -f 16 | tr '"' '\t'| cut -f 2  > $output/genomes_link

	# #######We online need the genomic.fna.gz 

	file="$output/genomes_link"
	while IFS= read -r line
	do
	    wget "$line/*_genomic.fna.gz" 
	done < "$file" 

	gzip -d *fna.gz 
	mv *genomic.fna $input
	rm *fna.gz

	######WE Separated genomes and plasmids
	genome_folder=$output/genomes_Shewanella/genomes
	plamid_folder=$output/genomes_Shewanella/plasmids
	mkdir -p genome_folder plamid_folder

	strains=( `ls $input` )
	for strain in "${strains[@]}"
	do
		genome_parser.rb $strain $genome_folder/$strain $plasmid_folder/$strain'_plasmids.fna'
	done
}

main_folder=`pwd`
. ~soft_bio_267/initializes/init_autoflow 

#####################################################
#######  MAIN 
##################################################
########1- genome reassigment (16S RNA, complete genome comparison and genome anotation) ####################  

project_path=`pwd`
data_path=$project_path'/data'
script_path=$project_path'/script'
genome_analysis=$SCRATCH/genome_shewanella
export PATH=$script_path:$PATH

if [ "$1" == "down" ]; then
 genomes_database data_path
fi

if [ "$1" == "char" ]; then
	. ~soft_bio_267/initializes/init_autoflow  

	# rm -r $genome_analysis
	# mkdir -p $genome_analysis
	# mkdir -p $data_path/total_genomes
	# ln -s $data_path/genomes_problem/*fasta $data_path/total_genomes
	# ln -s $data_path/genomes_Shewanella/genomes/*genomic.fna $data_path/total_genomes

	# ls $data_path/genomes_Shewanella/genomes > data/genome_candidate_list

	# TODO: mover todo lo q se escriba fuera de autoflow a su correspondiente carpeta
	vars=`echo "
		\\$data_path=$data_path,
		\\$ab1_name=$data_path/ab1_name,
		\\$prokaryotes_genomes_csv=$data_path/prokaryotes.csv,
		\\$genome_list=$data_path/genomes_link,
		\\$scripts_path=$script_path,
		\\$genome_candidate_list=$data_path/genome_candidate_list,
		\\$out_file_sibelia=/mnt/scratch/users/pab_001_uma/oliastencio/reasig_sibelia_results,
		\\$output_pyani=/mnt/scratch/users/pab_001_uma/oliastencio,
		\\$out_put_annotation=/mnt/scratch/users/pab_001_uma/oliastencio/dfast_annotation
		" | tr -d [:space:]`

	AutoFlow -c 1 -s -w script/shewanella_genome_analysis.af -V $vars -t "1-10:00:00" -o $genome_analysis
fi

#######Transposons ################

#transposons_identification            

#######4-Genes present in pathogenic strains###########
# genes_identification

#######5-Genomic island ##########
#genomic_island
