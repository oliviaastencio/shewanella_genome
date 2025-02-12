#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

data_path=$1
genome_analysis_path=$2
input=$3
output=$4

rm -r $output
mkdir -p $output
table=( 'interrupt_names' 'transposase_names' )

for tab in "${table[@]}"
	do 
	cat $input/*/lista_to_fasta.rb_0000/tp_case/*/$tab | sort -u > $output/'all_'$tab

	while read genome 
	do
		mkdir -p $output/$genome
		cds=`cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/$genome/CDSs`
		cat $input/$genome/lista_to_fasta.rb_0000/tp_case/*/$tab > $output/$genome/$tab
		cat $input/$genome/lista_to_fasta.rb_0000/tp_case/Total_tp | tr '\n' '\t' > $output/$genome/Total_tp_tab
		awk '{ print $1 "\t" $2 "\t" '"$cds"'}' $output/$genome/Total_tp_tab | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | cut -f 1,2 > $output/$genome/Tp_absolute
		awk '{ print $1 "\t" $2 "\t" '"$cds"'}' $output/$genome/Total_tp_tab | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | cut -f 1,4 > $output/$genome/Tp_relative

    	while read ID 
    	do
    	cero=0
    	uno=1
    	name=`echo $ID`
    	count=`grep -c "$ID" $output/$genome/$tab`
        
    	if [[ "$count" -eq 0 ]]; then 
            echo -e "$name\t$cero" >> $output/$genome/"tab_$tab"
        else         
            echo -e "$name\t$uno" >> $output/$genome/"tab_$tab"
        fi
        
        if [[ "$count" -gt 0 ]]; then 
            echo -e "$name\t$count" >> $output/$genome/"${tab}_total"
        else         
            echo -e "$name\t$cero" >> $output/$genome/"${tab}_total"
        fi
    	done < $output/'all_'$tab

	sed -i "1i category\t$genome" $output/$genome/"tab_$tab"
	sed -i "1i category\t$genome" $output/$genome/"${tab}_total"
	done < $data_path/all_genome_list

merge_tabular.rb $output/*/"tab_$tab" | sed s'/.fasta//g' > $output/"tab_$tab"

done 

cat $output/*/Tp_absolute > $output/Tp_absolute
cat $output/*/Tp_relative > $output/Tp_relative