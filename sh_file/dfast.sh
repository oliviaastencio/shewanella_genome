#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out


data_path=$1 
output=$2	#genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser

while read genome 
do 
	cds=`cat $output/"$genome"_CDSs`
 	awk '{ print $1 "\t" $2 "\t" '"$cds"'}' $output'/'$genome'_cog_table.txt' | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | grep -v "category" > $output'/'$genome'_cog_table_relative.txt'
 	cat $output'/'$genome'_cog_table_relative.txt' | cut -f 1,4 > $output/$genome'_cog_table'
	sed -i "1i category\t$genome" $output/$genome'_cog_table'
done < $data_path/genome_name

while read genome 
   do 
    accession=`grep "$genome" $data_path/total_genomes/classes.txt | cut -f 2`
    name=`grep "$genome" $data_path/total_genomes/classes.txt | cut -f 3`
    cds=`cat $output/"$accession".fna_CDSs`

  	awk '{ print $1 "\t" $2 "\t" '"$cds"'}' $output/$accession'.fna_cog_table.txt' | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | grep -v "category"> $output/$accession'.fna_cog_table_relative.txt'
  	cat $output/$accession'.fna_cog_table_relative.txt' | cut -f 1,4 > $output/$accession'.fna_cog_table'
  	sed -i "1i category\t$name" $output/$accession'.fna_cog_table'
done < $data_path/RefSeq_Accession.tsv
merge_tabular.rb $output/*cog_table > $output/Total_table.txt
grep "category" $output/Total_table.txt > $output/category_name
grep -v "category" $output/Total_table.txt > $output/COG_annotation_value
merge_tabular.rb  $data_path/COG_categories_all $output/COG_annotation_value | cut -f 2-133 > $output/COG_annotation_complete
cat $output/category_name $output/COG_annotation_complete | sed s'/.fasta//g' | sed s'/e_Shewanella_putrefaciens_//g' | sed s'/e_//g' | sed s'/_micro12//g' | sed s'/_micro13//g' | sed s'/_micro9//g' | sed s'/_micro22//g' | sed s'/_micro1//g' | sed s'/_1//g' > $output/COG_annotation
