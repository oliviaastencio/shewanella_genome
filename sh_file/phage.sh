#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=26G
#SBATCH --time=5-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

# ##################################
# ## Prophage by PHASTEST
# ##################################
data_path=$1
genome_analysis_path=$2

module load virsorter/2.2.4
rm -r $genome_analysis_path/prophage
mkdir -p $genome_analysis_path/prophage

	while read genome 
	do 
		length=`cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/$genome/length`

		mkdir -p $genome_analysis_path/prophage/$genome
		virsorter run -w  $genome_analysis_path/prophage/$genome -i $data_path/total_genomes/$genome  --min-length 1500 -j 4 all
		grep -v -c "seqname" $genome_analysis_path/prophage/$genome/final-viral-score.tsv >  $genome_analysis_path/prophage/$genome/final_number
		sed -i "1i $genome" $genome_analysis_path/prophage/$genome/final_number
		cat $genome_analysis_path/prophage/$genome/final_number | tr '\n' '\t' > $genome_analysis_path/prophage/$genome/final_number_tab
		awk '{ print $1 "\t" $2 "\t" '"$length"'}' $genome_analysis_path/prophage/$genome/final_number_tab | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | cut -f 1,2 | sed s'/.fasta//g' > $genome_analysis_path/prophage/$genome/final_number_absolute
		awk '{ print $1 "\t" $2 "\t" '"$length"'}' $genome_analysis_path/prophage/$genome/final_number_tab | awk '{ print $1 "\t" $2 "\t" $3 "\t" $2*100/$3}' | cut -f 1,4 | sed s'/.fasta//g' > $genome_analysis_path/prophage/$genome/final_number_relative
	
	done < $data_path/all_genome_list
cat $genome_analysis_path/prophage/*/*final_number_absolute > $genome_analysis_path/prophage/Total_phage
cat $genome_analysis_path/prophage/*/*final_number_relative > $genome_analysis_path/prophage/Total_phage_relative

