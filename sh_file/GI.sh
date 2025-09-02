#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

data_path=$2
genome_analysis_path=$3
input=$genome_analysis_path/genomic_island/Island_Viewer4_results/file_scv 
output=$genome_analysis_path/genomic_island/genomic_island_results
	
rm -r $output
mkdir -p $output

###########################################################################
#############################Total GI analized and upload 
##########################################################################

if [ "$1"  == "UP" ]; then 
	HTTP_API_token="b344dbd6-9cf9-720b-c8b4-7b6723104ba8" 
	input=$genome_analysis_path/char/dfast_0000/genome_annotation
	output=$genome_analysis_path/genomic_island

	mkdir -p $output
	rm -r $output/Island_Viewer4_results
	mkdir -p $output/Island_Viewer4_results

	while read ref_genome
	do 
		genome=`echo $ref_genome | tr ' ' '\t'| cut -f 1` 
		ref_num=`echo $ref_genome | tr ' ' '\t'| cut -f 2` 
		curl -X POST -H '$HTTP_API_token' -Fref_accnum="$ref_num" -Fgenome_file=@$input/$genome'_file'/genome.gbk -Fgenome_name="$genome"_genome -Femail_addr="oliastenciogomez@gmail.com" -Fformat_type="GENBANK" https://www.pathogenomics.sfu.ca/islandviewer/rest/submit/ > $output/Island_Viewer4_results/"$genome"_submit
		grep "token" $output/Island_Viewer4_results/"$genome"_submit | sed s'/    "token": "//g'| sed s'/"//g' > $output/Island_Viewer4_results/"$genome"_jobtoken
	done < $data_path/ref_genome_GI
fi 


###########################################################################
#############################Total GI results download  
##########################################################################

if [ "$1" == "result" ]; then 

##Query the status of already submitted genomes. Job_token is a job ID, the number where upload de file. I'ts not the same to my token ID or HTTP API Token
##curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/job_token/ -H 'X-authtoken:your_authentication_token'
	output=$genome_analysis_path/genomic_island/Island_Viewer4_results

	rm -r $output/file_scv
	mkdir -p $output/file_scv

	while read genome 
	do
		jobtoken=`head -n1 $output/"$genome"_jobtoken`
		job_token=`echo $jobtoken`
		curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/$jobtoken/download/csv/ -H '$HTTP_API_token' > $output/"$genome".csv
		mv $output/"$genome".csv $output/file_scv
	done < $data_path/all_genome_list
fi 


###########################################################################
#############################Total GI clean and enrichment
##########################################################################

if [ "$1" == "clean" ]; then

strains_island=( `ls $input` )

for strain in "${strains_island[@]}"
do
	rm -r $output/$strain
	mkdir -p $output/$strain
	mkdir -p $output/$strain/GI

	island_filtre.rb $input/$strain $output/$strain/$strain'_length'
	grep 'Predicted by at least' $output/$strain/$strain'_length'|cut -f 1,2,3,6,7,11,12 > $output/$strain/$strain'_Integrated'
	sed -i "1i Island_start\tIsland_end\tLength\tGene_ID\tLocus\tProduct\tExternal_Annotations" $output/$strain/$strain'_Integrated'
	grep -v "Island_start" $output/$strain/$strain'_Integrated' | cut -f 1,2,5 | awk '{ print "GI_"$1"_"$2 "\t" $3}' > $output/$strain/$strain'_LOCUS'
	echo $strain > $output/$strain/Total_genomic_island
	grep -v 'Island_start' $output/$strain/$strain'_Integrated' | cut -f1 | sort -u | wc -l >> $output/$strain/Total_genomic_island
	cat $output/$strain/$strain'_LOCUS' | cut -f 2 > $output/$strain/$strain'_locus'
	cat $output/$strain/$strain'_LOCUS' | cut -f 1 | sort -u > $output/$strain/$strain'_GI'
done 
 	
#######################################################################
 ### enrichment genomic island
########################################################################

while read genome 
	do 
		length=`cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/$genome/length`
		cat $output/$genome'.csv'/Total_genomic_island | tr '\n' '\t' | cut -f 1,2 >> $output/$genome'.csv'/Total_genomic_island_tab
		awk '{ print $1 "\t" $2 "\t" '"$length"'}' $output/$genome'.csv'/Total_genomic_island_tab | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | cut -f 1,4 >> $output/$genome'.csv'/Total_genomic_island_relative
    
	############################
	##create the  LOCUS list

	while read line 
	do	
		locus=`cut -f 1,3 $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/$genome/dfast_parser.txt | sed s'/LOCUS/_/g' | grep "$line" | cut -f 1 | sed s'/__/LOCUS_/g'`
		echo $line'tb'$locus | tr 'tb' '\t' | tr ' ' '\t' >> $output/$genome.csv/$genome.csv'_LOCUS_all'  
	done < $data_path/COG_categories

	############################
	### create the GI-LOCUS list
	while read GI 
	do 
		grep "$GI" $output/$genome.csv/$genome.csv'_LOCUS' | cut -f 2 | paste -s -d , >> $output/$genome.csv/GI/$GI'_GI_LOCUS'
		sed -i "1i $GI" $output/$genome.csv/GI/$GI'_GI_LOCUS'
		cat $output/$genome.csv/GI/$GI'_GI_LOCUS' | tr '\n' '\t' > $output/$genome.csv/GI/$GI'_GI_LOCUS_list'
		cat $output/$genome.csv/GI/*_GI_LOCUS_list | sed s'/GI/\nGI/g' > $output/$genome.csv/$genome'.csv_LOCUS_list'
		grep "GI" $output/$genome.csv/$genome.csv'_LOCUS_list' >  $output/$genome.csv/$genome'.csv_GI_LOCUS_list' 
	
	 done < $output/$genome.csv/$genome.csv'_GI'

	############################
	### enrichment cluster 

	 cd $output/$genome.csv/
	 . ~soft_bio_267/initializes/init_degenes_hunter
	 clusters_to_enrichment.R -i $genome'.csv_GI_LOCUS_list' --custom $genome.csv'_LOCUS_all' --funsys "" --gmt_id "" -k ""
	 cd ../

	if [[ -f "$output/$genome.csv/functional_enrichment/enr_cls_$genome.csv_LOCUS_all.csv" ]]; then
    awk '{if ($10<0.05) {print $0}}' "$output/$genome.csv/functional_enrichment/enr_cls_$genome.csv_LOCUS_all.csv" > "$output/$genome.csv/enr_cls_$genome.csv_LOCUS_all"
	else
    	> "$output/$genome.csv/enr_cls_$genome.csv_LOCUS_all"
	fi

	total_GI=`grep -c "GI_" $output/$genome.csv/enr_cls_$genome'.csv_LOCUS_all'`
	
	############################################
	### Tablet of enrichment and COG categories
	while read line 
	do 
	 	number=`grep -v "ID" $output/$genome.csv/enr_cls_$genome'.csv_LOCUS_all' | cut -f 2 | grep -c "$line"`
	 	echo $line't'$number | tr 't' '\t' >> $output/$genome.csv/'enrichment_GI_'$genome
	done < $data_path/COG_categories

	sed -i "1i category\t$genome" $output/$genome.csv/'enrichment_GI_'$genome
	############################################
	### Relative Tablet of enrichment and COG categories

	grep -v "category" "$output/$genome.csv/enrichment_GI_$genome" | awk -v total="$total_GI" '{perc=(total==0?0:$2*100/total); print $1 "\t" $2 "\t" perc}' | cut -f 1,3 > "$output/$genome.csv/relative_enrichment_GI_$genome"
	sed -i "1i category\t$genome" $output/$genome.csv/'relative_enrichment_GI_'$genome

	############################################
	## rm not necesary file
	rm -r $output/$genome.csv/GI $output/$genome.csv/$genome.csv* $output/$genome.csv/Total_genomic_island $output/$genome.csv/*LOCUS* $output/$genome.csv/*locus*

done < $data_path/all_genome_list


###########################################
## general TAB of all genomes

merge_tabular.rb $output/*/enrichment_GI_* | sed s'/.fasta//g' > $output/enrichment_GI_category
merge_tabular.rb $output/*/relative_enrichment_GI_* | sed s'/.fasta//g' > $output/enrichment_GI_category_relative

cat $output/*/*Total_genomic_island_tab | sed s'/.fasta.csv//g' > $output/../Total_GI
cat $output/*/*Total_genomic_island_relative | sed s'/.fasta.csv//g' > $output/../Total_GI_relative

fi