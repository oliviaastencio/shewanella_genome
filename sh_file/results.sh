#! /usr/bin/env bash

data_path=$1
results_path=$2
genome_analysis_path=$3

rm -r $results_path
mkdir -p $results_path

############################################
######### 16SRNA BLAST results #############
############################################

grep -h -v '#' $genome_analysis_path/char/blastn_0000/blast_16S/blast_* | awk '{if (($3>98.7)&&($15=100)) {print $0}}'| tr ' ' '\t' | sed s'/1_S/\tS/g' | cut -f 1,2,3,15 | sed s'/.fasta//g'| sed s'/_1492R//g' | sed s'/_/ /g'> $results_path/blast_16

############################################
######### Complete genome comparison #######
############################################

grep -v "Q_Shewanella_Pdp11" $genome_analysis_path/char/pyani_0000/genome_pyani_anim/matrix_identity_1.tab | grep -v "R_scaffold" | cut -f 1-9 | sed s'/Shewanella / S./g' | sort >> $results_path/pyani_identity
sed -i '1ishewanella strains \t Pdp11 \t SH12 \t SH16 \t SH4 \t SH6 \t SH9 \t SdM1 \t SdM2' $results_path/pyani_identity
grep -v "Q_Shewanella_Pdp11" $genome_analysis_path/char/pyani_0000/genome_pyani_anim/matrix_coverage_1.tab | grep -v "R_scaffold" | cut -f 1-9 | sed s'/Shewanella / S./g' | sort >> $results_path/pyani_coverage
sed -i '1ishewanella strains \t Pdp11 \t SH12 \t SH16 \t SH4 \t SH6 \t SH9 \t SdM1 \t SdM2' $results_path/pyani_coverage 
grep -v "shewanella strains" $results_path/pyani_identity | sed s'/ /_/g' | awk '{if (($2>0.9)||($3>0.9)||($4>0.9)||($5>0.9)||($6>0.9)||($7>0.9)||($8>0.9)||($9>0.9)) {print $0}}' | sed s'/_/ /g' | cut -f 1 > $results_path/pyani_identity_top

############################################
######### Synteny block and Sibelia #######
############################################

cp $genome_analysis_path/char/Sibelia_0000/Shewanella_Pdp11.fasta/Shewanella_baltica_128.fasta/circos/circos.png $results_path/S_baltica_128:Pdp11.png 
cp $genome_analysis_path/char/Sibelia_0000/Shewanella_Pdp11.fasta/Shewanella_putrefaciens_4H.fasta/circos/circos.png $results_path/S_putrefaciens_4H:Pdp11.png


############################################
######### Synteny block and Sibelia #######
############################################

### Genome annotation by DFAST
cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/*/COG_percent_tab |  sed s'/\tS/\nS/g' | sed s'/\tun/\nun/g' | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/COG_percent
sed -i '1ishewanella strains \t %' $results_path/COG_percent	
grep "category" $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/Total_cog_table | sed s'/Shewanella_/S./g'| sed s'/_/ /g' | tr '\t' '\n' > $results_path/name_list
grep "category" $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/Total_cog_table | sed s'/Shewanella_/S./g'| sed s'/_/ /g' > $results_path/name_tab
grep -v "category" $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/Total_cog_table | sed s'/Shewanella_/S./g'| sed s'/_/ /g' > $results_path/COG_annotation_data

while read line
do 
	if [ "$line" == "S.Pdp11.fasta" ]; then
		echo 1 | tr '\n' '\t' >> $results_path/number  
	fi 
	if [ "$line" != "S.Pdp11.fasta" ]; then
		echo 0 | tr '\n' '\t'>> $results_path/number  
	fi 
done < $results_path/name_list
cat $results_path/name_tab $results_path/number $results_path/COG_annotation_data | sed s'/RNA/\nRNA/g' > $results_path/COG_annotation
cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/Total_cog_table_relative | sed s'/Shewanella_/S./g'| sed s'/_/ /g' > $results_path/COG_annotation_relative
rm $results_path/name_list $results_path/number $results_path/COG_annotation_data

############################################
######### Transposon #######################
############################################

awk '{if ($2>0) {print $0}}' $genome_analysis_path/transposon/executions/Total_absolute | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/Total_absolute
sed -i '1ishewanella strains \t number' $results_path/Total_absolute  
awk '{if ($2>0) {print $0}}' $genome_analysis_path/transposon/executions/Total_relative | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/Total_relative
sed -i '1ishewanella strains \t %' $results_path/Total_relative

less -S $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,2,3,4 | tr ',' '\t' | cut -f 1,2,3,4 > $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/Pdp11_tp_1
less -S $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,5 | tr ',' '\t' | cut -f 1,2 > $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/Pdp11_tp_2
merge_tabular.rb $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/Pdp11_tp_1 $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/Pdp11_tp_2 | sed s'/Q_//g' | sed s'/_/ /' > $results_path/Pdp11_tp

############################################
######### Transposon matrix  ###############
############################################

cat $genome_analysis_path/transposon/executions/Tab_interrupt | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/Tp_interrupt
cat $genome_analysis_path/transposon/executions/Tab_interrupt_number | sed s'/Shewanella_/S./g' | sed s'/_/ /g'>  $results_path/Tp_interrupt_total
cat $genome_analysis_path/transposon/executions/Tab_transposase | sed s'/Shewanella_/S./g' | sed s'/_/ /g'>  $results_path/Tp_transposable
cat $genome_analysis_path/transposon/executions/Tab_transposase_number  | sed s'/Shewanella_/S./g' | sed s'/_/ /g'> $results_path/Tp_transposable_total

############################################
######### Tarsynflow #######################
############################################
 
grep -v 'Entry' $genome_analysis_path/genes_identification/Tarsynflow/results/Ref_specific_filtered_annotated_prots | cut -f 1,5,7 > $results_path/specific_genes


############################################
######### Genomic Island Total #############
############################################

cat $genome_analysis_path/genomic_island/genomic_island_results/Total_GI | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/GI_total
sed -i '1ishewanella strains \t GIs number' $results_path/GI_total

############################################
######### Genomic Island COG categories ####
############################################

grep "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_name
grep -v "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_value
merge_tabular.rb  $data_path/COG_categories_all $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_value | cut -f 2-139 > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI
cat $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_name $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/enrichment_GI_cat

############################################
######### Genomic Island COG categories ####
############################################

grep "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_name
grep -v "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_value
merge_tabular.rb  $data_path/COG_categories_all $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_value | cut -f 2-139 > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_relative
cat $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_name $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_relative | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/enrichment_GI_cat_relative


############################################
######### Prophage #########################
############################################

#cat $genome_analysis_path/prophage/Total_phage_top| sed s'/Shewanella/S./g' > $results_path/Total_phage
#sed -i '1ishewanella strains \t prophages number' $results_path/Total_phage
