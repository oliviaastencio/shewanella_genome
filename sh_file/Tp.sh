#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out


#### data 
data_path=$1
genome_analysis_path=$2
output=$genome_analysis_path/transposon/executions   #$genome_analysis_path/transposon/executions


##########################################################
################# TP matrix 
#########################################################
##interrupted genes
while read genome
	do 
		rm $output/$genome/all_interrupt_names $output/$genome/all_interrupt_names_tab
		while read ID 
		do 
			cero=0
			uno=1
			name=`echo $ID`
			strain=`grep $ID $output/$genome/lista_to_fasta.rb_0000/tp_case/all_interrupt_names`
			if [ "$name" == "$strain" ];then 
				echo "$name,$uno" >> $output/$genome/all_interrupt_names
			else 		 
				echo "$name,$cero" >> $output/$genome/all_interrupt_names
			fi
		done < $output/all_interrupt_names
	 	cat $output/$genome/all_interrupt_names | tr ',' '\t' > $output/$genome/all_interrupt_names_tab
done < $data_path/all_genome_list

##transposable elements 
while read genome 
	do
		rm $output/$genome/all_transposase_names $output/$genome/all_transposase_names_tab 
		while read ID 
		do 
			cero=0
			uno=1
			name=`echo $ID`
			strain=`grep $ID $output/$genome/lista_to_fasta.rb_0000/tp_case/all_transposase_names`
			if [ "$name" == "$strain" ];then 
				echo "$name,$uno" >> $output/$genome/all_transposase_names 
			fi 
			if [ "$name" != "$strain" ];then  
				echo "$name,$cero" >> $output/$genome/all_transposase_names 
			fi
		done < $output/all_transposase_names
		cat $output/$genome/all_transposase_names | tr ',' '\t' > $output/$genome/all_transposase_names_tab 	
done < $data_path/all_genome_list


##########################################################
################# TP file re-name 
#########################################################

rm $genome_analysis_path/transposon/executions/Total_tp
while read genome 
do 
	cp $output/$genome/lista_to_fasta.rb_0000/tp_case/Total_tp $output/$genome/lista_to_fasta.rb_0000/tp_case/Total_tp_clean
	sed -i "1i $genome" $output/$genome/lista_to_fasta.rb_0000/tp_case/Total_tp_clean
	cat $output/$genome/lista_to_fasta.rb_0000/tp_case/Total_tp_clean | tr '\n' '\t' | cut -f 1,2 > $output/$genome/lista_to_fasta.rb_0000/tp_case/Total_tp_tab
	sed -i "1i ID\t$genome" $output/$genome/all_interrupt_names_tab
	sed -i "1i ID\t$genome" $output/$genome/all_transposase_names_tab
done < $data_path/genome_name

while read genome 
   do 
    accession=`grep "$genome" $data_path/total_genomes/classes.txt | cut -f 2`
    name=`grep "$genome" $data_path/total_genomes/classes.txt | cut -f 3`
    cp $output/$accession.fna/lista_to_fasta.rb_0000/tp_case/Total_tp $output/$accession.fna/lista_to_fasta.rb_0000/tp_case/Total_tp_clean
  	sed -i "1i $name" $output/$accession.fna/lista_to_fasta.rb_0000/tp_case/Total_tp_clean
 	cat $output/$accession.fna/lista_to_fasta.rb_0000/tp_case/Total_tp_clean | tr '\n' '\t' | cut -f 1,2 > $genome_analysis_path/transposon/executions/$accession.fna/lista_to_fasta.rb_0000/tp_case/Total_tp_tab
    sed -i "1i ID\t$name" $output/$accession.fna/all_interrupt_names_tab
	sed -i "1i ID\t$name" $output/$accession.fna/all_transposase_names_tab
done < $data_path/RefSeq_Accession.tsv


cat $output/*/lista_to_fasta.rb_0000/tp_case/Total_tp_tab | sed s'/.fasta//g' | sed s'/e_Shewanella_putrefaciens_//g' | sed s'/e_//g' | sed s'/_micro12//g' | sed s'/_micro13//g' | sed s'/_micro9//g' | sed s'/_micro22//g' | sed s'/_micro1//g' | sed s'/_1//g' > $genome_analysis_path/transposon/executions/Total_tp
cat $genome_analysis_path/transposon/executions/Total_tp | sed s'/ /_/g' | awk '{if ($2>0) {print $0}}' | sed s'/_/ /g' > $genome_analysis_path/transposon/executions/Total_tp_top
merge_tabular.rb $output/*/all_interrupt_names_tab | sed s'/.fasta//g' | sed s'/e_Shewanella_putrefaciens_//g' | sed s'/e_//g' | sed s'/_micro12//g' | sed s'/_micro13//g' | sed s'/_micro9//g' | sed s'/_micro22//g' | sed s'/_micro1//g' | sed s'/_1//g' > $output/all_interrupt_names_tab
merge_tabular.rb $output/*/all_transposase_names_tab | sed s'/.fasta//g' | sed s'/e_Shewanella_putrefaciens_//g' | sed s'/e_//g' | sed s'/_micro12//g' | sed s'/_micro13//g' | sed s'/_micro9//g' | sed s'/_micro22//g' | sed s'/_micro1//g' | sed s'/_1//g' > $output/all_transposase_names_tab
 