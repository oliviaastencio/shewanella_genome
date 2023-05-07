#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
main_folder=`pwd`
. ~soft_bio_267/initializes/init_autoflow 

###############################################
##########transposon comparative##########
#########################################

function transposon_comparative {

	transposons_identification
	PATH=~soft_bio_267/programs/x86_64/scripts/:$PATH
	summary=/mnt/scratch/users/pab_001_uma/oliastencio/transposon_result/executions/e_Pdp11_1/transposons_finder.rb_0000/results/summary.txt  #shewanella_genomes_paper/data_script/transposons/tp_finder_results/results/summary.txt
	protein_fasta=data/tp_data/total_prots.fasta   #shewanella_genomes_paper/data_script/transposons/data/total_prots.fasta
	out_put=/mnt/scratch/users/pab_001_uma/oliastencio/transposon_result   #shewanella_genomes_paper/result/transposon/comparative

	#transposon complete sequence = disrupted protein fragment + transposon sequence + disrupted protein fragment
	mkdir -p $out_put/comparative
	module load blast_plus/2.2.30+
	module load cdhit/4.5.4
	awk '{ print $2":"$3 }' $summary > $out_put/comparative/tranposase_coordinates_length

	main_folder=`pwd`

	while read sequence
	do 
		folder='tp_'`echo $sequence | tr ':' '_'`
		call_folder=`echo $folder`
		mkdir -p $out_put/comparative/$folder

			while read coordinate
			do 
				echo '>'$folder'_'$coordinate > $out_put/comparative/$folder/transposase_"$coordinate".fasta
				length=`echo $coordinate | tr ':' '-'`
				head -n2 $out_put/tp_case/e_Pdp11_1/$call_folder/tp.fasta | cut -c $length | awk 'NF' >> $out_put/comparative/$folder/transposase_"$coordinate".fasta
			done < $out_put/comparative/tranposase_coordinates_length
			cat $out_put/comparative/$folder/transposase_* > $out_put/comparative/$folder/transposase.fasta
			grep $sequence $summary | cut -f 4 | tr ',' '\n' | tr '\t' '\n' > $out_put/comparative/$folder/protein_code
			lista_to_fasta.rb $protein_fasta $out_put/comparative/$folder/protein_code > $out_put/comparative/$folder/protein.fasta
			#We have used the transposon proteins sequence. 
			cd-hit -i $out_put/comparative/$folder/protein.fasta -o $out_put/comparative/$folder/protein_new.fasta -c 0.95

	done <$out_put/tp_case/e_Pdp11_1/tp_transposons_finder_coordinates

		mkdir -p $out_put/comparative/blastx_comparative
		####comparativa de las 3 parejas de tranposones similares por coordenadas.... tp_850951_856278 y tp_850951_856278; tp_1906106_1911827 y tp_1906106_1911827; también incluimos una tercera pareja aunque no tienen resultados similares 
		blastn -query $out_put/comparative/tp_850951_856278/transposase.fasta -subject $out_put/comparative/tp_1261471_1266798/transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_1261471_1266798
		awk '{if ($3>=35) {print $0}}' $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_1261471_1266798 | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary_tp_850951_856278_tp_1261471_1266798
		blastn -query $out_put/comparative/tp_1906106_1911827/transposase.fasta -subject $out_put/comparative/tp_3073685_3079406/transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_tp_1906106_1911827_tp_3073685_3079406
		awk '{if ($3>=35) {print $0}}' $out_put/comparative/blastx_comparative/blast_tp_1906106_1911827_tp_3073685_3079406 | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary_tp_1906106_1911827_tp_3073685_3079406
		blastn -query $out_put/comparative/tp_2693976_2700579/transposase.fasta -subject $out_put/comparative/tp_4349020_4355949/transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_tp_2693976_2700579_tp_4349020_4355949
		awk '{if ($3>=35) {print $0}}' $out_put/comparative/blastx_comparative/blast_tp_2693976_2700579_tp_4349020_4355949 | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary_tp_2693976_2700579_tp_4349020_4355949

		#####comparativa entre diferentes grupos
		blastn -query $out_put/comparative/tp_850951_856278/transposase.fasta -subject $out_put/comparative/tp_1906106_1911827/transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_1906106_1911827
		awk '{if ($3>=35) {print $0}}' $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_1906106_1911827 | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary_tp_850951_856278_tp_1906106_1911827
		blastn -query $out_put/comparative/tp_850951_856278/transposase.fasta  -subject $out_put/comparative/tp_2693976_2700579/transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_2693976_2700579
		awk '{if ($3>=35) {print $0}}' $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_2693976_2700579 | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary_tp_850951_856278_tp_2693976_2700579
}


##################################################################
##################Download Island_Viewer4_results
#################################################################

function Island_Viewer4_results {

	##Query the status of already submitted genomes. Job_token is a job ID, the number where upload de file. I'ts not the same to my token ID or HTTP API Token
	##curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/job_token/ -H 'X-authtoken:your_authentication_token'

	output=/mnt/scratch/users/pab_001_uma/oliastencio/genomic_island
	mkdir -p $output/Island_Viewer4_results/file_scv

	genomes=( 'e_Pdp11_1' 'e_Shewanella_putrefaciens_SH12_micro1' 'e_Shewanella_putrefaciens_SH4_micro9' 'e_Shewanella_putrefaciens_SH6_micro12' 'e_Shewanella_putrefaciens_SH9_micro13' 'e_Shewanella_putrefaciens_SH16_micro22' 'e_Shewanella_putrefaciens_SdM1' 'e_Shewanella_putrefaciens_SdM2' ) 
	HTTP_API_token=82950949-c400-0f86-dff6-07bf9adde409
	for genome in "${genomes[@]}"
	do 
	jobtoken=`head -n1 $output/Island_Viewer4_results/"$genome"_jobtoken`
	job_token=`echo $jobtoken`
	curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/$jobtoken/download/csv/ -H '$HTTP_API_token' > $output/Island_Viewer4_results/"$genome".csv
	mv $output/Island_Viewer4_results/"$genome".csv $output/Island_Viewer4_results/file_scv
	done
}


#####################################################################
##########################GENOMIC ISLAND FILTRE
####################################################################

function genomic_island_filter {
	
	PATH=~/micro_lab3/Pdp11_genome_analysis/genome_shewanella_Pdp11/script:$PATH
	export PATH
	module load ruby/2.4.1
	input=$output/Island_Viewer4_results/file_scv  
	output=/mnt/scratch/users/pab_001_uma/oliastencio/genomic_island 
	main_folder=`pwd`

	mkdir -p $output/genomic_island_results
	strains_island=( `ls $input` )
	for strain in "${strains_island[@]}"
	do
		mkdir -p $output/genomic_island_results/$strain
		island_filtre.rb $input/$strain $output/genomic_island_results/$strain/$strain'_length'
		grep 'Predicted by at least' $output/genomic_island_results/$strain/$strain'_length'|cut -f 1,2,3,6,7,11,12 > $output/genomic_island_results/$strain/$strain'_Integrated'
		sed -i "1i Island_start\tIsland_end\tLength\tGene_ID\tLocus\tProduct\tExternal_Annotations" $output/genomic_island_results/$strain/$strain'_Integrated'
		cat $output/genomic_island_results/$strain/$strain'_Integrated' | cut -f1 | sort -u | wc -l  > $output/genomic_island_results/$strain/Total_genomic_island
	done 
	cd $main_folder
}


#####################################################
#######  MAIN 
##################################################

#transposon_comparative          
Island_Viewer4_results   
genomic_island_filter

