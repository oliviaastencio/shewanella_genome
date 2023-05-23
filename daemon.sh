#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

######################################
#####FUNCTIONS
#####################################

##################################################
######Shewanella complete genomes database
#################################################
function genomes_database {
	output=$1
	module load ruby/2.4.1 
	rm -r $output/genomes_down
	mkdir -p $output/genomes_down/sequence_download
	
	main_folder=`pwd` 

	prokaryotes_genomes_csv=$output/prokaryotes.csv  #all shewanella genomes from https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/shewanella
	input=$output/genomes_down/sequence_download #sequence_download
	# rm $output/genomes_link
	cat $prokaryotes_genomes_csv | tr ',' '\t'| cut -f 16 | tr '"' '\t'| cut -f 2  > $output/genomes_link

	# #######We online need the genomic.fna.gz 

	file="$output/genomes_link"
	while IFS= read -r line
	do
	    wget "$line/*_genomic.fna.gz" 
	done < "$file" 

	gzip -d *fna.gz 
	mv *genomic.fna $input
	rm *fna.gz $input/*cds_from_genomic.fna $input/*rna_from_genomic.fna
	
	######WE Separated genomes and plasmids
	genome_folder=$output/genomes_down/genomes
	plasmid_folder=$output/genomes_down/plasmids
	mkdir -p $genome_folder $plasmid_folder

	strains=( `ls $input` )
	for strain in "${strains[@]}"
	do
		genome_parser.rb $input/$strain $genome_folder/$strain $plasmid_folder/$strain
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
genome_analysis_path=$SCRATCH/genome_analysis
transposon_analysis_path=$SCRATCH/transposon_result
export PATH=$script_path:$PATH
export PATH=$script_path'/scripts_transposon:'$PATH
export PATH=$script_path'/scripts_Tarsynflow:'$PATH

if [ "$1" == "down" ]; then  #####LISTO##############

	genomes_database
fi
 	
if [ "$1" == "ab1_clean" ]; then   #####LISTO##############
	ab1_clean.sh
fi

if [ "$1" == "char" ]; then  #####TODOS LOS BLOQUES LISTO########
	. ~soft_bio_267/initializes/init_autoflow  

	rm -r $genome_analysis_path/pyani_0000 $data_path/total_genomes $data_path/all_genome_list
	mkdir -p $data_path/total_genomes
	
	ln -s $data_path/genomes_problem/*fasta $data_path/total_genomes
	ln -s $data_path/genomes_down/genomes/*genomic.fna $data_path/total_genomes

	ls $data_path/genomes_down/genomes > $data_path/genome_candidate_list

	ls $data_path/genomes_problem > $data_path/all_genome_list
	ls $data_path/genomes_down/genomes >> $data_path/all_genome_list
	
	vars=`echo "
		\\$data_path=$project_path'/data',
		\\$ab1_name=$data_path/ab1_name,
		\\$genome_list=$data_path/genomes_link,
		\\$scripts_path=$script_path
		" | tr -d [:space:]`

	AutoFlow -c 1 -s -w $script_path/shewanella_genome_analysis.af -V $vars -t "1-10:00:00" -o $genome_analysis_path
fi

if [ "$1" == "protein_db" ]; then  #####LISTO##############

	protein_shewanella_db.sh 
fi 

if [ "$1" == "isescan" ]; then  #####LISTO##############

	. ~pedro/software/ISEScan/initialize 

	rm -r prediction $genome_analysis_path/transposon/executions/proteome $genome_analysis_path/transposon/executions/hmm

	nohup cat $genome_analysis_path/transposon/executions/fna.list | xargs -n 1 -P 8 -I{} isescan.py {} $genome_analysis_path/transposon/executions/proteome $genome_analysis_path/transposon/executions/hmm

fi 

if [ "$1" == "tp_case" ]; then  #####LISTO##############REVISAR RESULTADOS	
	
	while read line 
	do 
	mkdir -p $genome_analysis_path/transposon/executions/$line
	. ~soft_bio_267/initializes/init_autoflow 
	vars=`echo "
		\\$isscan_coordinates=$genome_analysis_path/transposon/executions/$line/ISSCAN_genome_coordinates,
		\\$genome_seq=$data_path/genomes_problem/$line.fasta,
		\\$prot_database=$data_path/tp_data
		" | tr -d [:space:]`
		AutoFlow -c 1 -s -w $script_path/tpflow -V $vars -o "$genome_analysis_path/transposon/executions/"$line $2

	done < $data_path/genome_name
fi 

if [ "$1" == "tp_comparative" ]; then   #####pendiente de ejecutar ####OPCIONAL!!

	tp_comparative.sh 

fi 

if [ "$1" == "genes" ]; then

	mkdir -p $genome_analysis_path/genes_identification/Tarsynflow/proteome

	reduce_prot_redundancy.sh

fi

if [ "$1" == "genes_comps" ]; then

	make_all_comps.sh

fi

if [ "$1" == "genes_results" ]; then
	get_all_results.sh
	extract_seqs.sh
fi 

 
if [ "$1"  == "GI" ]; then 

	genomic_island.sh
	
fi 

if [ "$1" == "GI_result" ]; then 

##Query the status of already submitted genomes. Job_token is a job ID, the number where upload de file. I'ts not the same to my token ID or HTTP API Token
##curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/job_token/ -H 'X-authtoken:your_authentication_token'

	mkdir -p $genome_analysis_path/genomic_island/Island_Viewer4_results/file_scv

	while read genome 
	do
		HTTP_API_token=82950949-c400-0f86-dff6-07bf9adde409
		jobtoken=`head -n1 $genome_analysis_path/genomic_island/Island_Viewer4_results/"$genome"_jobtoken`
		job_token=`echo $jobtoken`
		curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/$jobtoken/download/csv/ -H '$HTTP_API_token' > $genome_analysis_path/genomic_island/Island_Viewer4_results/"$genome".csv
		mv $genome_analysis_path/genomic_island/Island_Viewer4_results/"$genome".csv $genome_analysis_path/genomic_island/Island_Viewer4_results/file_scv
	done < $data_path/genome_name
fi 

if [ "$1" == "GI_filtre" ]; then 

	module load ruby/2.4.1
	input=$genome_analysis_path/genomic_island/Island_Viewer4_results/file_scv 
	output=$genome_analysis_path/genomic_island/genomic_island_results

	mkdir -p $output
	strains_island=( `ls $input` )
	for strain in "${strains_island[@]}"
	do
		mkdir -p $output/$strain
		island_filtre.rb $input/$strain $output/$strain/$strain'_length'
		grep 'Predicted by at least' $output/$strain/$strain'_length'|cut -f 1,2,3,6,7,11,12 > $output/$strain/$strain'_Integrated'
		sed -i "1i Island_start\tIsland_end\tLength\tGene_ID\tLocus\tProduct\tExternal_Annotations" $output/$strain/$strain'_Integrated'
		cat $output/$strain/$strain'_Integrated' | cut -f1 | sort -u | wc -l  > $output/$strain/Total_genomic_island
	done 
fi 



