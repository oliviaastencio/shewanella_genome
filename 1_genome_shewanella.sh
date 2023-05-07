#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
main_folder=`pwd`
. ~soft_bio_267/initializes/init_autoflow 
######################################
#####FUNCTIONS
#####################################

##################################################
######Shewanella complete genomes database 
#################################################
function genomes_database {
	16s_blastn_remote   
	rm -r data/genomes_Shewanella
	mkdir -p data/genomes_Shewanella
	mkdir -p data/genomes_Shewanella/sequence_download
		
	module load ruby/2.4.1 
	PATH=~oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/script:$PATH
	export PATH

	main_folder=`pwd` 

	 prokaryotes_genomes_csv=data/prokaryotes.csv  #all shewanella genomes from https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/shewanella
	 input=data/genomes_Shewanella/sequence_download #sequence_download
	 output=`pwd` #shewanella_genomes_paper/reassignment/genomes_Shewanella
	 rm data/genomes_link
	 cat $prokaryotes_genomes_csv | tr ',' '\t'| cut -f 16 | tr '"' '\t'| cut -f 2  > data/genomes_link

	# #######We online need the genomic.fna.gz 

	file="data/genomes_link"
	while IFS= read -r line
	do
	    wget "$line/*1_genomic.fna.gz" 
	    wget "$line/*2_genomic.fna.gz"
	    wget "$line/*3_genomic.fna.gz"
	    wget "$line/*4_genomic.fna.gz"
	done < "$file" 

	gzip -d *fna.gz 
	mv *genomic.fna $input
	rm *fna.gz

	######WE Separated genomes and plasmids

	mkdir -p data/genomes_Shewanella/genomes
	mkdir -p data/genomes_Shewanella/plasmids

	cd $input
	strains=( `ls` )
	for strain in "${strains[@]}"
	do
		genome_parser.rb $strain ~oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/data/genomes_Shewanella/genomes/$strain ~oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/data/genomes_Shewanella/plasmids/$strain'_plasmids.fna'
	done
	cd $main_folder
}
##################################################################################
##############16S_asignation###genome_asignation###########genome annotation#####
####################################################################################

function genomic_characteristics {
	. ~soft_bio_267/initializes/init_autoflow  

	output=/mnt/scratch/users/pab_001_uma/oliastencio/genome_shewanella
	# rm -r $output
	# mkdir -p $output
	# mkdir -p data/total_genomes
	# cp data/genomes_problem/*fasta data/total_genomes
	# cp data/genomes_Shewanella/genomes/*genomic.fna   data/total_genomes

	# ls /mnt/home/users/pab_001_uma/oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/data/genomes_Shewanella/genomes > data/genome_candidate_list

	vars=`echo "
				\\$data_path=/mnt/home/users/pab_001_uma/oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/data,
				\\$ab1_name=/mnt/home/users/pab_001_uma/oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/data/ab1_name,
				\\$prokaryotes_genomes_csv=/mnt/home/users/pab_001_uma/oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/data/prokaryotes.csv,
				\\$genome_list=/mnt/home/users/pab_001_uma/oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/data/genomes_link,
				\\$scripts_path=/mnt/home/users/pab_001_uma/oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/script,
				\\$genome_candidate_list=/mnt/home/users/pab_001_uma/oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/data/genome_candidate_list,
				\\$out_file_sibelia=/mnt/scratch/users/pab_001_uma/oliastencio/reasig_sibelia_results,
				\\$output_pyani=/mnt/scratch/users/pab_001_uma/oliastencio,
				\\$out_put_annotation=/mnt/scratch/users/pab_001_uma/oliastencio/dfast_annotation
				" | tr -d [:space:]`

	AutoFlow -c 1 -s -w script/shewanella_genome_analysis.af -V $vars -t "1-10:00:00" -o $output
}

##############################################
################TRANSPOSON_IDENTIFICATION#####
###############################################

function transposons_identification {

	framework_dir=`dirname $0`
	export CODE_PATH=$(readlink -f $framework_dir )
	export PATH=$CODE_PATH'/script/scripts_transposon:'$PATH
	 . ~soft_cvi_114/initializes/init_fln
	 
	taxonomy=22 #shewanella
	output=/mnt/scratch/users/pab_001_uma/oliastencio/transposon_result
	# rm -r $output
	# mkdir -p $output
	# mkdir -p $output/executions
	# mkdir -p data/tp_data/
	# export prot_fasta='data/tp_data/total_prots.fasta'

	curl  'https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=Shewanella' > $prot_fasta  #'https://www.uniprot.org:443/uniprot/?query=%20taxonomy:'$taxonomy'&format=fasta'
	cd data/tp_data
	pwd
	which make_user_db.rb
	make_user_db.rb -l -f total_prots.fasta -n local_database
	cd ../..
	rm $output/executions/fna.list
	rm -r prediction
	
	while read line ; do
		genome_seq=data/genomes_problem/$line.fasta
		echo $genome_seq >> $output/executions/fna.list
	done < data/genome_name
	
	#./script/isescan.sh
	
	. ~pedro/software/ISEScan/initialize 
	nohup cat $output/executions/fna.list | xargs -n 1 -P 8 -I{} isescan.py {} $output/executions/proteome $output/executions/hmm
 
	while read line 
	do 
	mkdir -p $output/executions/$line
	grep "insertion_sequence" prediction/genomes_problem/$line.fasta.gff | cut -f 1,4,5 > $output/executions/$line/ISSCAN_genome_coordinates 
	. ~soft_bio_267/initializes/init_autoflow 
	vars=`echo "
		\\$isscan_coordinates=$output/executions/$line/ISSCAN_genome_coordinates,
		\\$genome_seq=/mnt/home/users/pab_001_uma/oliastencio/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/data/genomes_problem/$line.fasta,
		\\$prot_database=$CODE_PATH/data/tp_data,
		\\$data=$CODE_PATH/data,
		\\$tp_case=/mnt/scratch/users/pab_001_uma/oliastencio/transposon_result/tp_case/$line
		" | tr -d [:space:]`
		AutoFlow -c 1 -s -w script/tpflow_1 -V $vars -o "$output/executions/"$line $2
	done < data/genome_name
}


####################################################
##########GENES SH IDENTIFICATION############
###############################################

function genes_identification {
	./script/get_seqs.sh
	./script/reduce_prot_redundancy.sh
	./script/make_all_comps.sh
	./script/get_all_results.sh
	./script/extract_seqs.sh
}

####################################################
#################GENOMIC ISLAND
####################################################

function genomic_island {

 	##1-login firts and I identify my HTTP API Token or my token ID and if this it not expired
	##2-upload de file where indicated the my token ID, email

 	#####################################complete genome######################################3
 	input=/mnt/scratch/users/pab_001_uma/oliastencio/dfast_annotation
 	output=/mnt/scratch/users/pab_001_uma/oliastencio/genomic_island
 	rm -r $output
	mkdir -p $output
	mkdir -p $output/Island_Viewer4_results
	HTTP_API_token=0597451f-c036-7302-c61d-2c66cd254e97

	genomes=( 'e_Pdp11_1' ) 
	
	for genome in "${genomes[@]}"
	do 
	curl -X POST -H '$HTTP_API_token' -Fgenome_file=@$input/$genome/genome.gbk -Fgenome_name=""$genome"_genome" -Femail_addr="oliastenciogomez@gmail.com" -Fformat_type="GENBANK" https://www.pathogenomics.sfu.ca/islandviewer/rest/submit/ > $output/Island_Viewer4_results/"$genome"_submit
	grep "token" $output/Island_Viewer4_results/"$genome"_submit | sed s'/    "token": "//g'| sed s'/"//g' > $output/Island_Viewer4_results/"$genome"_jobtoken
	done


	#######################################GENOME WITH Shewanella.sp FDAARGOS_354 of reference################################
	genomes=( 'e_Shewanella_putrefaciens_SH12_micro1' 'e_Shewanella_putrefaciens_SH4_micro9' )
	
	for genome in "${genomes[@]}"
	do 
	curl -X POST -H '$HTTP_API_token' -Fref_accnum="NZ_CP022089.2" -Fgenome_file=@$input/$genome/genome.gbk -Fgenome_name=""$genome"_genome" -Femail_addr="oliastenciogomez@gmail.com" -Fformat_type="GENBANK" https://www.pathogenomics.sfu.ca/islandviewer/rest/submit/ > $output/Island_Viewer4_results/"$genome"_submit
	grep "token" $output/Island_Viewer4_results/"$genome"_submit | sed s'/    "token": "//g'| sed s'/"//g' > $output/Island_Viewer4_results/"$genome"_jobtoken
	done
	#######################################GENOME WITH Shewanella sp. WE21 of reference###########################

	genomes=( 'e_Shewanella_putrefaciens_SH6_micro12' 'e_Shewanella_putrefaciens_SH9_micro13' 'e_Shewanella_putrefaciens_SH16_micro22' 'e_Shewanella_putrefaciens_SdM1' )

	for genome in "${genomes[@]}"
	do 
	curl -X POST -H '$HTTP_API_token' -Fref_accnum="NZ_CP023019.1" -Fgenome_file=@$input/$genome/genome.gbk -Fgenome_name=""$genome"_genome" -Femail_addr="oliastenciogomez@gmail.com" -Fformat_type="GENBANK" https://www.pathogenomics.sfu.ca/islandviewer/rest/submit/ > $output/Island_Viewer4_results/"$genome"_submit
	grep "token" $output/Island_Viewer4_results/"$genome"_submit | sed s'/    "token": "//g'| sed s'/"//g' > $output/Island_Viewer4_results/"$genome"_jobtoken
	done

	##############################GENOME with Shewanella putrefaciens WS13 of reference#################

	genomes=( 'e_Shewanella_putrefaciens_SdM2' )

	for genome in "${genomes[@]}"
	do 
	curl -X POST -H '$HTTP_API_token' -Fref_accnum="NZ_CP028435.1" -Fgenome_file=@$input/$genome/genome.gbk -Fgenome_name=""$genome"_genome" -Femail_addr="oliastenciogomez@gmail.com" -Fformat_type="GENBANK" https://www.pathogenomics.sfu.ca/islandviewer/rest/submit/ > $output/Island_Viewer4_results/"$genome"_submit
	grep "token" $output/Island_Viewer4_results/"$genome"_submit | sed s'/    "token": "//g'| sed s'/"//g' > $output/Island_Viewer4_results/"$genome"_jobtoken
	done
 }


#####################################################
#######  MAIN 
##################################################
########1- genome reassigment (16S RNA, complete genome comparison and genome anotation) ####################  

#genomes_database
genomic_characteristics

#######Transposons ################

#transposons_identification            

#######4-Genes present in pathogenic strains###########
# genes_identification

#######5-Genomic island ##########
#genomic_island


