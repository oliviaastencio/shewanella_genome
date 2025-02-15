#! /usr/bin/env bash

data_path=$1
results_path=$2
genome_analysis_path=$3

rm -r $results_path
mkdir -p $results_path
mkdir -p $results_path/report_img

########################################################
####data clean##########################################
genus=$4
initial=$(echo "$genus" | cut -f 1)
reference_genome=$5

linea_num=0

while read -r line; do
    name=$(echo "$line" | sed -e 's/.fasta//g' -e "s/${genus}_//g")
    ((linea_num++))  

    declare "strain$linea_num=$name"
done < "$data_path/genome_name"

echo "Our analyzed strains"
for ((i=1; i<=linea_num; i++)); do
    eval "echo strain$i=\$strain$i"
done

###### Define the sed command correctly in a variable
######renaming Shewanella unculture, naming all Shewanella strains by their initial "S"
sed_command_lowbar="sed -e 's/uncultured_${genus}_sp/${genus}_sp._uncultured/g' \
                -e 's/${genus}_/${initial}._/g' \
                -e 's/${initial}._sp._/${genus}_sp._/g'"
#####renaming our strains from "Shewanellas" to "Shewanella sp."
sed_name="sed -e 's/${genus}.sp./${genus}_sp./g' \
				-e 's/${initial}._${strain1}/${genus}_sp._${strain1}/g' \
				-e 's/${initial}._${strain2}/${genus}_sp._${strain2}/g' \
				-e 's/${initial}._${strain3}/${genus}_sp._${strain3}/g' \
				-e 's/${initial}._${strain4}/${genus}_sp._${strain4}/g' \
				-e 's/${initial}._${strain5}/${genus}_sp._${strain5}/g' \
				-e 's/${initial}._${strain6}/${genus}_sp._${strain6}/g' \
				-e 's/${initial}._${strain7}/${genus}_sp._${strain7}/g' \
				-e 's/${initial}._${strain8}/${genus}_sp._${strain8}/g'"

############################################
######### 16SRNA BLAST results #############
############################################

grep -h -v '#' $genome_analysis_path/char/blastn_0000/blast_16S/blast_* | tr ' ' '_' | awk '{if (($3>98.7)&&($15=100)) {print $0}}'| tr ' ' '\t' | cut -f 1,3,15,16 | sed s'/_1492R//g' | sed s'/.fasta//g' | tr '_' ' '  > $results_path/report_img/blast_16 

############################################
######### Complete genome comparison #######
############################################

##########################################
############ rename the results of pyani

while read label 
do 
	ID=`echo "$label" | cut -f 2`
	header=`echo "$label" | cut -f 3`
	table=( 'matrix_identity' 'matrix_coverage' ) 

	for tab in "${table[@]}"
	do
		tail -n +2 $genome_analysis_path/char/pyani_0000/genome_pyani_anim/$tab'_1.tab'| grep "$header" | sed s"/${header}/${ID}/g" | cut -f 1-9  >> $results_path/'pyani'_$tab
	done
done < $data_path/total_genomes/classes.txt


##############################################################################################
############ rename as known or unknown species in data of report

table=( 'matrix_identity' 'matrix_coverage' ) 

header="${genus} strains\t${initial}._${strain1}\t${initial}._${strain4}\t${initial}._${strain5}\t${initial}._${strain6}\t${initial}._${strain7}\t${initial}._${strain8}\t${initial}._${strain2}\t${initial}._${strain3}"

for tab in "${table[@]}"
do
	echo -e "$header" | cat - $results_path/'pyani'_$tab > temp && mv temp $results_path/'pyani'_$tab
	cat $results_path/'pyani'_$tab | eval "$sed_command_lowbar" | sort | grep -Ev "$strain1|$strain2|$strain3|$strain4|$strain5|$strain6|$strain7|$strain8" >> $results_path/report_img/'pyani'_$tab
	echo -e "$header" | cat - $results_path/report_img/'pyani'_$tab | eval "$sed_name" | sed s'/_/ /g' > temp && mv temp $results_path/report_img/'pyani'_$tab
done


########################################################################################################
######### Synteny block and Sibelia ####################################################################
########################################################################################################
############ take the genome with the highest identity by pyani and extract its png image from Sibelia

while read genome 
do
	name=`echo $genome | sed -e "s/${genus}_//g" -e s'/.fasta//g'`
	col=`head $results_path/report_img/pyani_matrix_identity | tr '\t' '\n' | nl | grep "^.*$name" | awk '{print $1}'`
	reference=`cat $results_path/pyani_matrix_identity | cut -f 1,$col | tr ':' '\t' | cut -f 1,3 | grep -Ev "${genus} strains|$strain1|$strain2|$strain3|$strain4|$strain5|$strain6|$strain7|$strain8" | sort -k 2 | tail -n 1 | cut -f 1 `
	cp $genome_analysis_path/char/Sibelia_0000/$genome/*$reference'.fasta'/circos/circos.png  $results_path/report_img/$name'_'$reference'.png'
	cp $genome_analysis_path/char/Sibelia_0000/$genome/$reference_genome'.fasta'/circos/circos.png  $results_path/report_img/$name'_'$reference_genome'.png'
done < $data_path/genome_name

############################################
######### Genome annotation by DFAST #######
############################################
	

table=( 'Total_cog_table' 'Total_cog_table_relative' ) 

for tab in "${table[@]}"
do
	path=$genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser
	cat $path/$tab | sed s'/.fasta//g' > $results_path/$tab
done

##################################################################################
######### MOBILOME: Genomic Island, Prophage and Transposon #######################
###################################################################################
cat $genome_analysis_path/transposon/results/Tp_absolute | sed s'/.fasta//g' > $genome_analysis_path/transposon/Tp_absolute
cat $genome_analysis_path/transposon/results/Tp_relative | sed s'/.fasta//g' > $genome_analysis_path/transposon/Tp_relative


table=( 'Total_phage' 'Total_phage_relative' 'Total_GI' 'Total_GI_relative' 'Tp_absolute' 'Tp_relative' ) 

for tab in "${table[@]}"
do
	awk '{if ($2>0) {print $0}}' $genome_analysis_path/*/$tab | sed s'/ /_/g'| eval "$sed_command_lowbar" > $results_path/$tab
	echo -e "${genus} strains \t number" | cat - $results_path/$tab >  temp && mv temp $results_path/$tab
	grep -v "${genus} strains" $results_path/$tab | cut -f 1  > $results_path/$tab'_name'

	while read line
	do
		if [ "$line" == "${initial}._${strain1}" ]; then
			echo "$line" Pdp | tr ' ' '\t' >> $results_path/number  
		fi 
		if [ "$line" == "${initial}._${strain7}" ] || [ "$line" == "${initial}._${strain8}" ] || [ "$line" == "${initial}._${strain4}" ] || [ "$line" == "${initial}._${strain5}" ] || [ "$line" == "${initial}._${strain6}" ]; then
			echo "$line" SH | tr ' ' '\t' >> $results_path/number  
		fi
		if [ "$line" == "${initial}._${strain2}" ] || [ "$line" == "${initial}._${strain3}" ] ; then
			echo "$line" SdM | tr ' ' '\t' >> $results_path/number  
		fi
		if [ "$line" != "${initial}._${strain1}" ] && [ "$line" != "${initial}._${strain2}" ] && [ "$line" != "${initial}._${strain3}" ] && [ "$line" != "${initial}._${strain4}" ] && [ "$line" != "${initial}._${strain5}" ] && [ "$line" != "${initial}._${strain6}" ] && [ "$line" != "${initial}._${strain7}" ] && [ "$line" != "${initial}._${strain8}" ] ; then
			echo "$line" "${genus}" | tr ' ' '\t' >> $results_path/number  
		fi 
	done < $results_path/$tab'_name'

	grep -v "number" $results_path/$tab > $results_path/$tab'_data' 
	merge_tabular.rb $results_path/number $results_path/$tab'_data' > $results_path/$tab'_all'
	cat  $results_path/$tab'_all' | grep -v "${genus} strains" | sort -k3,3g -n | eval "$sed_name" > $results_path/report_img/'report_'$tab
	echo -e "${genus} \t strains \t number" | cat - $results_path/report_img/'report_'$tab | sed s'/_/ /g' >  temp && mv temp $results_path/report_img/'report_'$tab
	rm $results_path/$tab'_name' $results_path/number $results_path/$tab'_data' $results_path/$tab'_all'
done


#################################################################
######### Transposon of ${strain1} #############################
#################################################################

less -S $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,2,3,4 | tr ',' '\t' | cut -f 1,2,3,4 > $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_1
less -S $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,5 | tr ',' '\t' | cut -f 1,2 > $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_2
merge_tabular.rb $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_1 $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_2 | sed s'/Q_//g' | sed s'/_/ /' > $results_path/report_img/Tp_"$strain1"

tp_parser.rb $genome_analysis_path/transposon/results/tab_interrupt_names $results_path/Tp_interrupt
tp_parser.rb $genome_analysis_path/transposon/results/tab_transposase_names $results_path/Tp_transposase


############################################
######### Tarsynflow #######################
############################################
 
grep -v 'Entry' $genome_analysis_path/genes_identification/Tarsynflow/results/Ref_specific_filtered_annotated_prots | cut -f 1,5,7 > $results_path/report_img/specific_genes


#################################################################
######### Genomic Island COG categories enrichment (absolute) ####
table=( 'enrichment_GI_category' 'enrichment_GI_category_relative' ) 

for tab in "${table[@]}"
do
	path=$genome_analysis_path/genomic_island/genomic_island_results

	grep "category" $path'/'$tab > $path/$tab'_name'
	grep -v "category" $path'/'$tab > $path/$tab'_value'
	merge_tabular.rb  $data_path/COG_categories_all $path'/'$tab'_value' | cut -f 2-139 > $path/enrichment_GI
	cat $path'/'$tab'_name' $path/enrichment_GI > $results_path/$tab
	rm $path/$tab'_name' $path/$tab'_value' $path/enrichment_GI
done 


#####################################################################################
#### COG annotation, GI enrichment, interrupted and transposed protein ##############
#####################################################################################

table=( 'Total_cog_table' 'Total_cog_table_relative' 'enrichment_GI_category' 'enrichment_GI_category_relative' 'Tp_interrupt' 'Tp_transposase' ) 

for tab in "${table[@]}"
do
	grep "category" $results_path/$tab | eval "$sed_command_lowbar" | cut -f 2-138 | tr '\t' '\n' > $results_path/name_list 
	grep "category" $results_path/$tab | eval "$sed_command_lowbar" | sed s'/.fasta//g' > $results_path/name_tab 
	first_data=`grep -v "category" $results_path/$tab | cut -f 1 | head -n 1`
 	grep -v "category" $results_path/$tab  > $results_path/temp_data 
	while read line
	do 
		if [ "$line" == "${initial}._${strain1}" ]; then
			echo "    Pdp" | tr '\n' '\t' >> $results_path/number  
		fi 
		if [ "$line" == "${initial}._${strain7}" ] || [ "$line" == "${initial}._${strain8}" ] || [ "$line" == "${initial}._${strain4}" ] || [ "$line" == "${initial}._${strain5}" ] || [ "$line" == "${initial}._${strain6}" ]; then
			echo "    SH"| tr '\n' '\t' >> $results_path/number  
		fi 
		if [ "$line" == "${initial}._${strain2}" ] || [ "$line" == "${initial}._${strain3}" ] ; then
			echo "    SdM" | tr '\n' '\t' >> $results_path/number  
		fi
		if [ "$line" != "${initial}._${strain1}" ] && [ "$line" != "${initial}._${strain2}" ] && [ "$line" != "${initial}._${strain3}" ] && [ "$line" != "${initial}._${strain4}" ] && [ "$line" != "${initial}._${strain5}" ] && [ "$line" != "${initial}._${strain6}" ] && [ "$line" != "${initial}._${strain7}" ] && [ "$line" != "${initial}._${strain8}" ] ; then
			echo "${genus}" | tr '\n' '\t'>> $results_path/number  
		fi 
	done < $results_path/name_list

sed -i '1istrains' $results_path/number 
cat $results_path/number | tr '\n' '\t' > $results_path/number_tab
cat $results_path/name_tab $results_path/number_tab $results_path/temp_data  | sed s"/${first_data}/\n${first_data}/g" | eval "$sed_name" | sed s'/_/ /g' > $results_path/report_img/'report_'$tab  
rm  $results_path/name_list $results_path/name_tab $results_path/temp_data $results_path/number $results_path/number_tab 
done
