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
##############################################
#### Clean data output from IslandViewer 4
##############################################

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
done 
    
#######################################################################
  ### select the GI and its LOCUS with COG annotation ,  Rename the file and genomes
########################################################################

#######1) our genomes 
while read genome 
	do 
	cds=`cat $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/"$genome"_CDSs`

	### select the LOCUS with COG annotation and create the COG list 
	while read locus
	do
		grep "$locus" $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/dfast_parser_$genome.txt | cut -f 3 | sed s'/,//g' >> $output/$genome.csv/$genome'.csv_categories'
		grep "$locus" $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/dfast_parser_$genome.txt | cut -f 1 >> $output/$genome.csv/$genome'.csv_LOCUS_COG'
		grep "$locus" $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/dfast_parser_$genome.txt >> $output/$genome.csv/$genome'.csv_COG_list'
		GI_parser.rb $output/$genome.csv/$genome.csv_categories $output/$genome.csv/$genome.csv'_categories'
		awk '{ print $1 "\t" $2 "\t" '"$cds"'}' $output/$genome.csv/category_table.txt | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | grep -v "categories" | cut -f 1,4 > $output/$genome.csv/category_table
		sed -i "1i categories\t$genome" $output/$genome.csv/category_table
	done < $output/$genome.csv/$genome.csv_locus
	cat $output/$genome.csv/Total_genomic_island | tr '\n' '\t' | cut -f 1,2 >> $output/$genome.csv/Total_genomic_island_tab

### create the genome-GI list
	while read line
	do 
		grep "$line" $output/$genome.csv/$genome.csv'_LOCUS' | cut -f 1  >> $output/$genome.csv/$genome.csv'_GI_all'
	done < $output/$genome.csv/$genome.csv'_LOCUS_COG'
	cat $output/$genome.csv/$genome.csv'_GI_all' | sort -u > $output/$genome.csv/$genome.csv'_GI'
	cat $output/$genome.csv/$genome.csv'_GI' | paste -s -d , > $output/$genome.csv/$genome.csv'_genome_GI'
	sed -i "1i $genome" $output/$genome.csv/$genome.csv'_genome_GI'
	cat $output/$genome.csv/$genome.csv'_genome_GI' | tr '\n' '\t' > $output/$genome.csv/$genome.csv'_GI_list'

	### create the GI-LOCUS list
	while read GI 
	do 
		grep "$GI" $output/$genome.csv/$genome.csv'_LOCUS' | cut -f 2 | paste -s -d , >> $output/$genome.csv/GI/$GI'_GI_LOCUS'
		grep "$GI" $output/$genome.csv/$genome.csv'_LOCUS' | cut -f 1,2 >> $output/$genome.csv/GI/$GI'_GI_LOCUS_tab'
		grep "$GI" $output/$genome.csv/$genome.csv'_LOCUS' | cut -f 2 >> $output/$genome.csv/GI/$GI'_LOCUS_tab'
		sed -i "1i $GI" $output/$genome.csv/GI/$GI'_GI_LOCUS'
		while read COG 
		do 
			grep "$COG" $output/$genome.csv/$genome'.csv_COG_list' >> $output/$genome.csv/GI/$GI'_categories_list'
			grep "$COG" $output/$genome.csv/$genome'.csv_COG_list' | cut -f 1 >> $output/$genome.csv/GI/$GI'_categories_COG_list'
			while read cat
			do
				grep $cat $output/$genome.csv/GI/$GI'_GI_LOCUS_tab' >> $output/$genome.csv/GI/$GI'_GI_LOCUS_COG'
			done < $output/$genome.csv/GI/$GI'_categories_COG_list'
		done < $output/$genome.csv/GI/$GI'_LOCUS_tab'
	awk '{ print $2 "\t" $1}' $output/$genome.csv/GI/$GI'_GI_LOCUS_COG' | sort -u > $output/$genome.csv/GI/$GI'_GI_LOCUS_COG_tab'
	merge_tabular.rb $output/$genome.csv/GI/$GI'_categories_list' $output/$genome.csv/GI/$GI'_GI_LOCUS_COG_tab' | awk '{ print $1"_"$2":"$3}' | paste -s -d , > $output/$genome.csv/GI/$GI'_LOCUS_GI_list'
	sed -i "1i $GI" $output/$genome.csv/GI/$GI'_LOCUS_GI_list' 
	cat $output/$genome.csv/GI/$GI'_LOCUS_GI_list'| tr '\n' '\t' > $output/$genome.csv/GI/$GI'_LOCUS_list'
	done < $output/$genome.csv/$genome.csv'_GI'
	cat $output/$genome.csv/GI/*_LOCUS_list | sed s'/GI/\nGI/g' > $output/$genome.csv/$genome.csv'_LOCUS_list'
done < $data_path/genome_name

#######2) NCBI genomes

# select the LOCUS with COG annotation 
while read genome 
do 
accession=`grep "$genome" $data_path/total_genomes/classes.txt | cut -f 2`
name=`grep "$genome" $data_path/total_genomes/classes.txt | cut -f 3`
cds=`cat $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/"$accession".fna_CDSs`
	while read locus
	do
		grep "$locus" $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/dfast_parser_$accession.fna.txt | cut -f 3 | sed s'/,//g' >> $output/$accession.fna.csv/$accession.fna.csv_categories
		grep "$locus" $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/dfast_parser_$accession.fna.txt | cut -f 1 >> $output/$accession.fna.csv/$accession.fna'.csv_LOCUS_COG'
		grep "$locus" $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/dfast_parser_$accession.fna.txt >> $output/$accession.fna.csv/$accession.fna'.csv_COG_list'
		GI_parser.rb $output/$accession.fna.csv/$accession.fna.csv_categories $output/$accession.fna.csv/$accession.fna.csv'_categories'
		awk '{ print $1 "\t" $2 "\t" '"$cds"'}' $output/$accession.fna.csv/category_table.txt | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | grep -v "categories" | cut -f 1,4 > $output/$accession.fna.csv/category_table
		sed -i "1i categories\t$name" $output/$accession.fna.csv/category_table
	done < $output/$accession.fna.csv/$accession.fna.csv_locus
	sed -i "1i $name" $output/$accession.fna.csv/Total_genomic_island 
	cat $output/$accession.fna.csv/Total_genomic_island | tr '\n' '\t'| cut -f 1,3 >> $output/$accession.fna.csv/Total_genomic_island_tab

	### create the GI list
	while read line
	do 
		grep $line $output/$accession.fna.csv/$accession.fna.csv'_LOCUS' | cut -f 1  >> $output/$accession.fna.csv/$accession.fna.csv'_GI_all'
	done < $output/$accession.fna.csv/$accession.fna.csv'_LOCUS_COG'
	cat $output/$accession.fna.csv/$accession.fna.csv'_GI_all' | sort -u > $output/$accession.fna.csv/$accession.fna.csv'_GI'
	cat $output/$accession.fna.csv/$accession.fna.csv'_GI' | paste -s -d , > $output/$accession.fna.csv/$accession.fna.csv'_genome_GI'
	sed -i "1i $name" $output/$accession.fna.csv/$accession.fna.csv'_genome_GI'
	cat $output/$accession.fna.csv/$accession.fna.csv'_genome_GI' | tr '\n' '\t' > $output/$accession.fna.csv/$accession.fna.csv'_GI_list'

	### create the LOCUS list
	while read GI 
	do 
		grep "$GI" $output/$accession.fna.csv/$accession.fna.csv'_LOCUS' | cut -f 2 | paste -s -d , >> $output/$accession.fna.csv/GI/$GI'_GI_LOCUS'
		grep "$GI" $output/$accession.fna.csv/$accession.fna.csv'_LOCUS' | cut -f 1,2 >> $output/$accession.fna.csv/GI/$GI'_GI_LOCUS_tab'
		grep "$GI" $output/$accession.fna.csv/$accession.fna.csv'_LOCUS' | cut -f 2 >> $output/$accession.fna.csv/GI/$GI'_LOCUS_tab'
		sed -i "1i $GI" $output/$accession.fna.csv/GI/$GI'_GI_LOCUS'
		while read COG 
		do 
			grep "$COG" $output/$accession.fna.csv/$accession.fna'.csv_COG_list' >> $output/$accession.fna.csv/GI/$GI'_categories_list'
			grep "$COG" $output/$accession.fna.csv/$accession.fna'.csv_COG_list' | cut -f 1 >> $output/$accession.fna.csv/GI/$GI'_categories_COG_list'
			while read cat
			do
				grep $cat $output/$accession.fna.csv/GI/$GI'_GI_LOCUS_tab' >> $output/$accession.fna.csv/GI/$GI'_GI_LOCUS_COG'
			done < $output/$accession.fna.csv/GI/$GI'_categories_COG_list'
		done < $output/$accession.fna.csv/GI/$GI'_LOCUS_tab'
		awk '{ print $2 "\t" $1}' $output/$accession.fna.csv/GI/$GI'_GI_LOCUS_COG' | sort -u > $output/$accession.fna.csv/GI/$GI'_GI_LOCUS_COG_tab'
		merge_tabular.rb $output/$accession.fna.csv/GI/$GI'_categories_list' $output/$accession.fna.csv/GI/$GI'_GI_LOCUS_COG_tab' | awk '{ print $1":"$2":"$3}' | paste -s -d , > $output/$accession.fna.csv/GI/$GI'_LOCUS_GI_list'
		sed -i "1i $GI" $output/$accession.fna.csv/GI/$GI'_LOCUS_GI_list' 
		cat $output/$accession.fna.csv/GI/$GI'_LOCUS_GI_list'| tr '\n' '\t' > $output/$accession.fna.csv/GI/$GI'_LOCUS_list'
	done < $output/$accession.fna.csv/$accession.fna.csv'_GI'
	cat $output/$accession.fna.csv/GI/*_LOCUS_list | sed s'/GI/\nGI/g' > $output/$accession.fna.csv/$accession.fna.csv'_LOCUS_list'
done < $data_path/RefSeq_Accession.tsv
merge_tabular.rb $output/*/*category_table | sed s'/.fasta//g' | sed s'/e_Shewanella_putrefaciens_//g' | sed s'/e_//g' | sed s'/_micro12//g' | sed s'/_micro13//g' | sed s'/_micro9//g' | sed s'/_micro22//g' | sed s'/_micro1//g' | sed s'/_1//g'> $output/Total_category
cat $output/*/*tab | sed s'/.fasta.csv//g' | sed s'/e_Shewanella_putrefaciens_//g' | sed s'/e_//g' | sed s'/_micro12//g' | sed s'/_micro13//g' | sed s'/_micro9//g' | sed s'/_micro22//g' | sed s'/_micro1//g' | sed s'/_1//g'> $output/Total_GI
cat $output/Total_GI | sed s'/ /_/g' | awk '{if ($2>10) {print $0}}' | sed s'/_/ /g' > $output/Total_GI_top
rm -r $output/*/GI 
rm -r $output/*/*category_table.txt $output/*/*categories $output/*/*genome_GI $output/*/*GI_all