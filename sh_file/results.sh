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
genus_low_case=`echo $genus| tr 'S' 's'`
initial=S
shewanella_type1=Shewanella_baltica_128
shewanella_type2=Shewanella_putrefaciens_4H

scaffold=R_scaffold 
strain1=Pdp11
strain2=SH12
strain3=SH16
strain4=SH4
strain5=SH6
strain6=SH9
strain7=SdM1
strain8=SdM2

other=MAG:

###### Define the sed command correctly in a variable
sed_command_space="sed -e 's/uncultured ${genus}.sp/    ${genus}.sp.uncultured/g' \
                -e 's/${genus} /    ${initial}./g' \
                -e 's/${initial}.sp./${genus}.sp./g'"
sed_command_lowbar="sed -e 's/uncultured_${genus}_sp/${genus}_sp.uncultured/g' \
                -e 's/${genus}_/${initial}./g' \
                -e 's/${initial}.sp./${genus}.sp./g'"
sed_name="sed -e 's/${genus}.sp./${genus} sp. /g' \
				-e 's/${initial}.${strain1}/${genus} sp.${strain1}/g' \
				-e 's/${initial}.${strain2}/${genus} sp.${strain2}/g' \
				-e 's/${initial}.${strain3}/${genus} sp.${strain3}/g' \
				-e 's/${initial}.${strain4}/${genus} sp.${strain4}/g' \
				-e 's/${initial}.${strain5}/${genus} sp.${strain5}/g' \
				-e 's/${initial}.${strain6}/${genus} sp.${strain6}/g' \
				-e 's/${initial}.${strain7}/${genus} sp.${strain7}/g' \
				-e 's/${initial}.${strain8}/${genus} sp.${strain8}/g'"

############################################
######### 16SRNA BLAST results #############
############################################

grep -h -v '#' $genome_analysis_path/char/blastn_0000/blast_16S/blast_* | awk '{if (($3>98.7)&&($15=100)) {print $0}}'| tr ' ' '\t' | sed s'/1_S/\tS/g' | cut -f 1,2,3,15 | sed s'/.fasta//g'| sed s'/_1492R//g' | sed s'/_/ /g'> $results_path/report_img/blast_16

############################################
######### Complete genome comparison #######
############################################	
table=( 'matrix_identity' 'matrix_coverage' ) 

for tab in "${table[@]}"
do
	grep -v "Q_${genus}_${strain1}" $genome_analysis_path/char/pyani_0000/genome_pyani_anim/$tab'_1.tab' | grep -v "$scaffold" | cut -f 1-9 | sed s"/${other}/     /g" | eval "$sed_command_space" | sort >> $results_path/report_img/'pyani'_$tab 
	echo -e "${genus} strains\t${initial}.${strain1}\t${initial}.${strain2}\t${initial}.${strain3}\t${initial}.${strain4}\t${initial}.${strain5}\t${initial}.${strain6}\t${initial}.${strain7}\t${initial}.${strain8}" | cat - $results_path/report_img/'pyani'_$tab | eval "$sed_name" > temp && mv temp $results_path/report_img/'pyani'_$tab 
done 
############################################
######### Synteny block and Sibelia #######
############################################

cp $genome_analysis_path/char/Sibelia_0000/"$genus"_"$strain1".fasta/"$shewanella_type1".fasta/circos/circos.png $results_path/report_img/"$shewanella_type1".png 
cp $genome_analysis_path/char/Sibelia_0000/"$genus"_"$strain1".fasta/"$shewanella_type2".fasta/circos/circos.png $results_path/report_img/"$shewanella_type2".png


############################################
######### Synteny block and Sibelia #######
############################################

### Genome annotation by DFAST
cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/*/COG_percent_tab |  sed s'/\tS/\nS/g' | sed s'/\tun/\nun/g' | eval "$sed_command_lowbar"> $results_path/COG_percent
echo -e "${genus} strains \t %" | cat - $results_path/COG_percent | eval "$sed_name" >  temp && mv temp $results_path/COG_percent	
	

cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/Total_cog_table | sed s'/.fasta//g' > $results_path/Cog_table
cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/Total_cog_table_relative | sed s'/.fasta//g' > $results_path/Cog_table_relative


#################################################################
######### Transposon absolute and relative #######################
#################################################################

awk '{if ($2>0) {print $0}}' $genome_analysis_path/transposon/executions/Total_absolute | sed s'/ /_/g'| eval "$sed_command_lowbar" > $results_path/Total_tp_absolute
echo -e "${genus_low_case} strains \t number" | cat - $results_path/Total_tp_absolute >  temp && mv temp $results_path/Total_tp_absolute
awk '{if ($2>0) {print $0}}' $genome_analysis_path/transposon/executions/Total_relative | sed s'/ /_/g'| eval "$sed_command_lowbar" > $results_path/Total_tp_relative
echo -e "${genus_low_case} strains \t %" | cat - $results_path/Total_tp_relative >  temp && mv temp $results_path/Total_tp_relative

#################################################################
######### Transposon of ${strain1} (Pdp11)##################################
#################################################################

less -S $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,2,3,4 | tr ',' '\t' | cut -f 1,2,3,4 > $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_1
less -S $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,5 | tr ',' '\t' | cut -f 1,2 > $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_2
merge_tabular.rb $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_1 $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_2 | sed s'/Q_//g' | sed s'/_/ /' > $results_path/report_img/Tp_"$strain1"
cp $genome_analysis_path/transposon/executions/Tab_interrupt  $results_path/Tp_interrupt
cp $genome_analysis_path/transposon/executions/Tab_transposase $results_path/Tp_transposase

############################################
######### Tarsynflow #######################
############################################
 
grep -v 'Entry' $genome_analysis_path/genes_identification/Tarsynflow/results/Ref_specific_filtered_annotated_prots | cut -f 1,5,7 > $results_path/report_img/specific_genes


############################################
######### Genomic Island Total #############
############################################

cat $genome_analysis_path/genomic_island/genomic_island_results/Total_GI | sed s'/ /_/g'| eval "$sed_command_lowbar" > $results_path/Total_GI
echo -e "${genus_low_case} strains \t GIs number" | cat - $results_path/Total_GI >  temp && mv temp $results_path/Total_GI

#################################################################
######### Genomic Island COG categories enrichment (absolute) ####


grep "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_name
grep -v "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_value
merge_tabular.rb  $data_path/COG_categories_all $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_value | cut -f 2-139 > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI
cat $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_name $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI > $results_path/enrichment_GI


##################################################################
######### Genomic Island COG categories enrichment (relative) ####

grep "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_name
grep -v "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_value
merge_tabular.rb  $data_path/COG_categories_all $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_value | cut -f 2-139 > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_relative
cat $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_name $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_relative > $results_path/enrichment_GI_relative

############################################
######### Prophage #########################
############################################

cat Total_phage | sed s"/${strain1}/${initial}.${strain1}/g" | sed s"/${strain2}/${initial}.${strain2}/g" | sed s"/${strain3}/${initial}.${strain3}/g" | sed s"/${strain4}/${initial}.${strain4}/g" | sed s"/${strain5}/${initial}.${strain5}/g" | sed s"/${strain6}/${initial}.${strain6}/g" | sed s"/${strain7}/${initial}.${strain7}/g"  | sed s"/${strain8}/${initial}.${strain8}/g"  | sed s"/${initial}. /${initial}./g" | sed s'/ /_/g' | sed s"/${other}_//g" | sed s'/ /_/g'| eval "$sed_command_lowbar" | sed s"/${genus_low_case}_strains_/${genus_low_case} strains/g" > $results_path/Total_phage #| sed s'/_number/number/g'

##################################################################################################
#### COG annotation, GI enrichment, interrupted and transposed protein ##############
##################################################################################################

table=( 'Cog_table' 'Cog_table_relative' 'enrichment_GI' 'enrichment_GI_relative' 'Tp_interrupt' 'Tp_transposase' ) 

for tab in "${table[@]}"
do
	grep "category" $results_path/$tab | sed s'/_/ /g' | eval "$sed_command_space" | cut -f 2-138 | tr '\t' '\n' > $results_path/name_list 
	grep "category" $results_path/$tab | sed s'/_/ /g' | eval "$sed_command_space" | sed s'/.fasta//g' > $results_path/name_tab 
	grep -v "category" $results_path/$tab | sed s'/_/ /g' | eval "$sed_command_space"  > $results_path/temp_data 
	while read line
	do 
		if [ "$line" == "${initial}.${strain1}" ]; then
			echo Pdp | tr '\n' '\t' >> $results_path/number  
		fi 
		if [ "$line" == "${initial}.${strain2}" ] || [ "$line" == "${initial}.${strain3}" ] || [ "$line" == "${initial}.${strain4}" ] || [ "$line" == "${initial}.${strain5}" ] || [ "$line" == "${initial}.${strain6}" ]; then
			echo SH | tr '\n' '\t' >> $results_path/number  
		fi
		if [ "$line" == "${initial}.${strain7}" ] || [ "$line" == "${initial}.${strain8}" ] ; then
			echo SdM | tr '\n' '\t' >> $results_path/number  
		fi
		if [ "$line" != "${initial}.${strain1}" ] && [ "$line" != "${initial}.${strain2}" ] && [ "$line" != "${initial}.${strain3}" ] && [ "$line" != "${initial}.${strain4}" ] && [ "$line" != "${initial}.${strain5}" ] && [ "$line" != "${initial}.${strain6}" ] && [ "$line" != "${initial}.${strain7}" ] && [ "$line" != "${initial}.${strain8}" ] ; then
			echo "${genus}" | tr '\n' '\t'>> $results_path/number  
		fi 
	done < $results_path/name_list

sed -i '1istrains' $results_path/number 
cat $results_path/number | tr '\n' '\t' > $results_path/number_tab
cat $results_path/name_tab $results_path/number_tab $results_path/temp_data | sed s'/RNA/\nRNA/g' | sed s'/A0A6G9QKY9/\nA0A6G9QKY9/g' | sed s'/A0A1S6HLN8/\nA0A1S6HLN8/g' | eval "$sed_name" > $results_path/report_img/'report_'$tab  
rm  $results_path/name_list $results_path/name_tab $results_path/temp_data $results_path/number $results_path/number_tab  #$results_path/$tab
done

###########################################################################################################
#### Genomic Island, Prophage , Transposon absolute and Transposon relative  barplot with var_attr ######## 
##########################################################################################################


table=( 'Total_GI' 'Total_phage' 'Total_tp_absolute' 'Total_tp_relative' ) 

for tab in "${table[@]}"
do

	grep -v "${genus_low_case} strains" $results_path/$tab | cut -f 1  > $results_path/$tab'_name'

	while read line
	do
		if [ "$line" == "${initial}.${strain1}" ]; then
			echo "$line" Pdp | tr ' ' '\t' >> $results_path/number  
		fi 
		if [ "$line" == "${initial}.${strain2}" ] || [ "$line" == "${initial}.${strain3}" ] || [ "$line" == "${initial}.${strain4}" ] || [ "$line" == "${initial}.${strain5}" ] || [ "$line" == "${initial}.${strain6}" ]; then
			echo "$line" SH | tr ' ' '\t' >> $results_path/number  
		fi
		if [ "$line" == "${initial}.${strain7}" ] || [ "$line" == "${initial}.${strain8}" ] ; then
			echo "$line" SdM | tr ' ' '\t' >> $results_path/number  
		fi
		if [ "$line" != "${initial}.${strain1}" ] && [ "$line" != "${initial}.${strain2}" ] && [ "$line" != "${initial}.${strain3}" ] && [ "$line" != "${initial}.${strain4}" ] && [ "$line" != "${initial}.${strain5}" ] && [ "$line" != "${initial}.${strain6}" ] && [ "$line" != "${initial}.${strain7}" ] && [ "$line" != "${initial}.${strain8}" ] ; then
			echo "$line" "${genus}" | tr ' ' '\t' >> $results_path/number  
		fi 
	done < $results_path/$tab'_name'

grep -v "number" $results_path/$tab > $results_path/$tab'_data' 
merge_tabular.rb $results_path/number $results_path/$tab'_data' > $results_path/$tab'_all'
cat  $results_path/$tab'_all' | grep -v "${genus_low_case} strains" | sort -k3,3 -n | sed s'/_/ /g' > $results_path/report_img/'report_'$tab
echo -e "${genus} \t strains \t number" | cat - $results_path/report_img/'report_'$tab | eval "$sed_name" >  temp && mv temp $results_path/report_img/'report_'$tab
rm $results_path/$tab'_name' $results_path/number $results_path/$tab'_data' $results_path/$tab'_all'
done