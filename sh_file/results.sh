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
initial=$(echo "$genus" | cut -c1)
reference_S_putrefaciens=$5
reference_S_baltica=$6

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
strains=("${strain1}" "${strain2}" "${strain3}" "${strain4}" "${strain5}" "${strain6}" "${strain7}" "${strain8}")

sed_command_all="sed -e 's/uncultured_${genus}_sp/${genus}_sp._uncultured/g' -e 's/${genus}.sp./${genus}_sp./g'"

for s in "${strains[@]}"; do
  sed_command_all+=" -e 's/${genus}_${s}/${genus}_sp._${s}/g'"
done

# ## modify names##########################################
# ##For each strain, we keep the genus abbreviation (e.g., S.), the first four letters of the species name, and the strain identifier, omitting any duplicate strain IDs. This format makes visualization easier.
# ##For example, Shewanella baltica OS185 becomes S. balt OS185.

awk_cmd='awk -F"\t" '\''{
  n=split($1,p,"_"); seg=substr(p[2],1,4); resto="";
  for(i=3;i<=(n>=4?3:n);i++){ sub(/:.*/,"",p[i]); resto=resto?resto" "p[i]:p[i] }
  $1=substr(p[1],1,1)". "seg" "resto;
  printf "%s", $1; for(i=2;i<=NF;i++) printf "\t%s",$i; print "";
}'\'''

############################################
######### 16SRNA BLAST results #############
############################################

grep -h -v '#' $genome_analysis_path/char/blastn_0000/blast_16S/blast_* | tr ' ' '_' | awk '{if (($3>98.7)&&($15=100)) {print $0}}'| tr ' ' '\t' | cut -f 1,3,15,16 | sed s'/_1492R//g' | sed s'/.fasta//g' | tr '_' ' '  > $results_path/report_img/blast_16 

############################################
######### Complete genome comparison #######
############################################

##########################################
############ rename the results of pyani

tables=(matrix_identity matrix_coverage)
header="${genus} strains\t${initial}._${strain1}\t${initial}._${strain4}\t${initial}._${strain5}\t${initial}._${strain6}\t${initial}._${strain7}\t${initial}._${strain8}\t${initial}._${strain2}\t${initial}._${strain3}"

for tab in "${tables[@]}"; do
  pyani="$results_path/pyani_$tab"; report="$results_path/report_img/pyani_$tab"; > "$pyani"
  while IFS=$'\t' read -r _ ID hdr; do
    tail -n +2 "$genome_analysis_path/char/pyani_0000/genome_pyani_anim/${tab}_1.tab" | grep -F "$hdr" | sed "s/$hdr/$ID/g" | cut -f 1-9 >> "$pyani"   
  done < "$data_path/total_genomes/classes.txt"
  {  echo -e "$header"; cat "$pyani"; }  | sort | grep -Ev "$strain1|$strain2|$strain3|$strain4|$strain5|$strain6|$strain7|$strain8"  | eval "$sed_command_all" | eval "$awk_cmd" | sed "1i$header" | sed 's/_/ /g' > "$report"  
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
done < $data_path/genome_name
	cp $genome_analysis_path/char/Sibelia_0000/'Shewanella_'$strain1'.fasta'/$reference_S_putrefaciens'.fasta'/circos/circos.png  $results_path/report_img/$strain1'_'$reference_S_putrefaciens'.png'
	cp $genome_analysis_path/char/Sibelia_0000/'Shewanella_'$strain1'.fasta'/$reference_S_baltica'.fasta'/circos/circos.png  $results_path/report_img/$strain1'_'$reference_S_baltica'.png'


############################################
######### Genome annotation by DFAST #######
############################################

table=( 'Total_cog_table' 'Total_cog_table_relative' )
path="$genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser"

for tab in "${table[@]}"; do
  sed 's/.fasta//g' "$path/$tab" > "$results_path/$tab"
done

#################################################################
######### Transposon of ${strain1} #############################
#################################################################

less -S $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,2,3,4 | tr ',' '\t' | cut -f 1,2,3,4 > $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_1
less -S $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,5 | tr ',' '\t' | cut -f 1,2 > $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_2
merge_tabular.rb $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_1 $genome_analysis_path/transposon/executions/"$genus"_"$strain1".fasta/transposons_finder.rb_0000/results/"$strain1"_tp_2 | sed s'/Q_//g' | sed s'/_/ /' > $results_path/report_img/Tp_"$strain1"

tp_parser.rb $genome_analysis_path/transposon/results/tab_interrupt_names $results_path/Tp_interrupt
tp_parser.rb $genome_analysis_path/transposon/results/tab_transposase_names $results_path/Tp_transposase


#################################################################
######### Genomic Island COG categories enrichment (absolute) ####

table=( 'enrichment_GI_category' 'enrichment_GI_category_relative' )
path=$genome_analysis_path/genomic_island/genomic_island_results

for tab in "${table[@]}"; do
  grep "category" "$path/$tab" > "$path/${tab}_name"
  grep -v "category" "$path/$tab" > "$path/${tab}_value"
  merge_tabular.rb "$data_path/COG_categories_all" "$path/${tab}_value" | cut -f 2-139 > "$path/enrichment_GI" &&
  cat "$path/${tab}_name" "$path/enrichment_GI" > "$results_path/$tab" &&
  rm "$path/${tab}_name" "$path/${tab}_value" "$path/enrichment_GI"
done

######################################################################################################################################################
######### Heading level 2: Genomic Island, Prophage and Transposon, COG annotation, GI enrichment, interrupted and transposed protein #######################
######################################################################################################################################################


for f in Tp_absolute Tp_relative; do
  sed 's/.fasta//g' "$genome_analysis_path/transposon/results/$f" > "$genome_analysis_path/transposon/$f"
done

tables=(Total_phage Total_phage_relative Total_GI Total_GI_relative Tp_absolute Tp_relative Total_cog_table Total_cog_table_relative enrichment_GI_category enrichment_GI_category_relative Tp_interrupt Tp_transposase)
strain_pdp=("${genus}_sp._${strain1}")
strain_sh=("${genus}_sp._${strain7}" "${genus}_sp._${strain8}" "${genus}_sp._${strain4}" "${genus}_sp._${strain5}" "${genus}_sp._${strain6}")
strain_sdm=("${genus}_sp._${strain2}" "${genus}_sp._${strain3}")

assign_label() {
  [[ " ${strain_pdp[*]} " == *" $1 "* ]] && echo -ne "Pdp\t" && return
  [[ " ${strain_sh[*]} " == *" $1 "* ]] && echo -ne "SH\t" && return
  [[ " ${strain_sdm[*]} " == *" $1 "* ]] && echo -ne "SdM\t" && return
  echo -ne "${genus}\t"
}

for tab in "${tables[@]}"; do
  input="$results_path/$tab"
  if [[ $tab =~ ^(Total_phage|Total_phage_relative|Total_GI|Total_GI_relative|Tp_absolute|Tp_relative)$ ]]; then    ##### Total_phage Total_phage_relative Total_GI Total_GI_relative Tp_absolute Tp_relative
    awk '$2>0' "$genome_analysis_path"/*/"$tab" | sed 's/ /_/g' | eval "$sed_command_all" > "$input"  
    { echo -e "${genus} strains\t number"; cat "$input"; } > temp && mv temp "$input"
    grep -v "${genus} strains" "$input" | cut -f 1 > "$results_path/${tab}_name"
    : > "$results_path/number"
    while read -r line; do
      case "$line" in
        "${genus}_sp._${strain1}") echo -e "$line\tPdp" ;;  
        "${genus}_sp._${strain2}"| "${genus}_sp._${strain3}") echo -e "$line\tSdM" ;;
        "${genus}_sp._${strain4}"| "${genus}_sp._${strain5}"| "${genus}_sp._${strain6}"| "${genus}_sp._${strain7}"| "${genus}_sp._${strain8}") echo -e "$line\tSH" ;;
        *) echo -e "$line\t${genus}" ;;
      esac >> "$results_path/number"
    done < "$results_path/${tab}_name"

    grep -v number "$input" > "$results_path/${tab}_data"
    merge_tabular.rb "$results_path/number" "$results_path/${tab}_data" > "$results_path/${tab}_all"
    cat "$results_path/${tab}_all" | grep -v "${genus} strains" | sed 's/ /_/g' | eval "$sed_command_all" | sort -k3,3g -n | eval "$awk_cmd"  > "$results_path/report_img/report_${tab}"  # 
    { echo -e "${genus}\t strains\t number"; cat "$results_path/report_img/report_${tab}"; } | sed 's/_/ /g' > temp && mv temp "$results_path/report_img/report_${tab}"  
    rm "$results_path/${tab}_name" "$results_path/number" "$results_path/${tab}_data" "$results_path/${tab}_all"

  else   ### tab Total_cog_table Total_cog_table_relative enrichment_GI_category enrichment_GI_category_relative Tp_interrupt Tp_transposase
    grep "category" "$input" | cut -f 2-138 | eval "$sed_command_all" | tr '\t' '\n' > "${results_path}/name_list" 
    cat "${results_path}/name_list" | eval "$awk_cmd" | tr '\n' '\t' > "${results_path}/name_tab"  
    sed -i -e '1s/^/strains\t/' -e '$a\' "${results_path}/name_tab"
    first_data=$(grep -v "category" "$input" | cut -f 1 | head -n 1)
    grep -v "category" "$input" > "${results_path}/temp_data"
    > "${results_path}/number"
    while read -r line; do assign_label "$line"; done < "${results_path}/name_list" > "${results_path}/number"
    sed -i -e '1s/^/strains\t/' "${results_path}/number"
    cat "${results_path}/name_tab" "${results_path}/number" "${results_path}/temp_data" | sed "s/${first_data}/\n${first_data}/" | sed 's/_/ /g'  > "${results_path}/report_img/report_${tab}"
    rm "${results_path}"/{name_list,name_tab,temp_data,number}
  fi
done


############################################
######### Tarsynflow #######################
############################################
 
grep -v 'Entry' $genome_analysis_path/genes_identification/Tarsynflow/results/Ref_specific_filtered_annotated_prots | cut -f 1,5,7 > $results_path/report_img/specific_genes

