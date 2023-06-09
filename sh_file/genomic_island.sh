#! /usr/bin/env bash

##1-login firts and I identify my HTTP API Token or my token ID and if this it not expired
##2-upload de file where indicated the my token ID, email

 #####################################complete genome######################################3
HTTP_API_token=$1
input=$2
output=$3
rm -r $output
mkdir -p $output/Island_Viewer4_results

# TODO: crear a mano archivo tabulado con cepa a estudiar y cepa de referencia en segunda columna, leer en bucle para lanzar consultas
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