#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#module load ruby/2.4.1
data_path=$1
genome_analysis_path=$2
input=$genome_analysis_path/genomic_island/Island_Viewer4_results/file_scv 
output=$genome_analysis_path/genomic_island/genomic_island_results
	
rm -r $output
mkdir -p $output

#############################################
### Clean data output from IslandViewer 4
#############################################

strains_island=( `ls $input` )

for strain in "${strains_island[@]}"
do
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
  ### select the GI and its LOCUS with COG annotation ,  Rename the file and genomes
########################################################################

while read genome 
	do 
	cds=`cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/"$genome"/CDSs`
	while read line 
	do	
		locus=`cut -f 1,3 $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/$genome/dfast_parser.txt | sed s'/LOCUS/_/g' | grep "$line" | cut -f 1 | sed s'/__/LOCUS_/g'`
		echo $line'tb'$locus | tr 'tb' '\t' | tr ' ' '\t' >> $output/$genome.csv/$genome.csv'_LOCUS_all'  
	done < $data_path/COG_categories
	#grep "COG" $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/dfast_parser_$genome.txt | cut -f 1,3 > $output/$genome.csv/$genome.csv'_LOCUS_all'
	
	# select the LOCUS with COG annotation at genome and create the COG list 

	while read locus
	do
		grep "$locus" $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/$genome/dfast_parser.txt | cut -f 3 | sed s'/,//g' >> $output/$genome.csv/$genome'.csv_categories'
		GI_parser.rb $output/$genome.csv/$genome.csv_categories $output/$genome.csv/$genome.csv'_categories'
		awk '{ print $1 "\t" $2 "\t" '"$cds"'}' $output/$genome.csv/category_table.txt | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | grep -v "categories" | cut -f 1,2 > $output/$genome.csv/category_table
		awk '{ print $1 "\t" $2 "\t" '"$cds"'}' $output/$genome.csv/category_table.txt | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | grep -v "categories" | cut -f 1,4 > $output/$genome.csv/category_table_relative
	done < $output/$genome.csv/$genome.csv_locus
	sed -i "1i categories\t$genome" $output/$genome.csv/category_table
	sed -i "1i categories\t$genome" $output/$genome.csv/category_table_relative
	cat $output/$genome.csv/Total_genomic_island | tr '\n' '\t' | cut -f 1,2 >> $output/$genome.csv/Total_genomic_island_tab
	### create the GI-LOCUS list
	while read GI 
	do 
		grep "$GI" $output/$genome.csv/$genome.csv'_LOCUS' | cut -f 2 | paste -s -d , >> $output/$genome.csv/GI/$GI'_GI_LOCUS'
		sed -i "1i $GI" $output/$genome.csv/GI/$GI'_GI_LOCUS'
		cat $output/$genome.csv/GI/$GI'_GI_LOCUS' | tr '\n' '\t' > $output/$genome.csv/GI/$GI'_GI_LOCUS_list'
		cat $output/$genome.csv/GI/*_GI_LOCUS_list | sed s'/GI/\nGI/g' > $output/$genome.csv/$genome'.csv_LOCUS_list'
		grep "GI" $output/$genome.csv/$genome.csv'_LOCUS_list' >  $output/$genome.csv/$genome'.csv_GI_LOCUS_list' 
	
	 done < $output/$genome.csv/$genome.csv'_GI'

	 cd $output/$genome.csv/
	 . ~soft_bio_267/initializes/init_degenes_hunter
	 clusters_to_enrichment.R -i $genome'.csv_GI_LOCUS_list' --custom $genome'.csv_LOCUS_all' --funsys "" --gmt_id "" -k ""
	 cd ../

	 awk '{if ($10<0.05) {print $0}}' $output/$genome.csv/functional_enrichment/enr_cls_$genome'.csv_LOCUS_all.csv' > $output/$genome.csv/enr_cls_$genome'.csv_LOCUS_all'
	 total_GI=`grep -c "GI_" $output/$genome.csv/enr_cls_$genome'.csv_LOCUS_all'`

	 while read line 
	 do 
	 	number=`grep -v "ID" $output/$genome.csv/enr_cls_$genome'.csv_LOCUS_all' | cut -f 2 | grep -c "$line"`
	 	echo $line't'$number | tr 't' '\t' >> $output/$genome.csv/'enrichment_GI_'$genome
	 done < $data_path/COG_categories
	 sed -i "1i category\t$genome" $output/$genome.csv/'enrichment_GI_'$genome
	 grep -v "category" $output/$genome.csv/'enrichment_GI_'$genome | awk '{ print $1 "\t" $2 "\t" '"$total_GI"'}' | awk '{ print $1 "\t" $2 "\t" $2*100/$3}' | cut -f 1,3 > $output/$genome.csv/'relative_enrichment_GI_'$genome
	 sed -i "1i category\t$genome" $output/$genome.csv/'relative_enrichment_GI_'$genome
done < $data_path/all_genome_list


#### Tab Genomic Island 
merge_tabular.rb $output/*/*category_table | sed s'/.fasta//g' > $output/Total_category
merge_tabular.rb $output/*/*category_table_relative | sed s'/.fasta//g' > $output/Total_category_relative

#### Tab Category table genomic Island
merge_tabular.rb $output/*/enrichment_GI_* | sed s'/.fasta//g' > $output/enrichment_GI_category
merge_tabular.rb $output/*/relative_enrichment_GI_* | sed s'/.fasta//g' > $output/enrichment_GI_category_relative

cat $output/*/*tab | sed s'/.fasta.csv//g' > $output/Total_GI
cat $output/Total_GI | sed s'/ /_/g' | awk '{if ($2>10) {print $0}}' | sed s'/_/ /g' > $output/Total_GI_top
rm -r $output/*/GI 
rm -r $output/*/*category_table.txt $output/*/*categories #$output/*/*GI_all $output/*/*genome_GI

#GI-LOCUS: LOCUS_list: _LOCUS_category
#LOCUS-CATEGORIES:  _GI_LOCUS_list