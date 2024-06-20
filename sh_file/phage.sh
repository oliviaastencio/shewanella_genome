#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

data_path=$1
output=$2  #genome_analysis_path/prophage
rm $output/*/phage_number_tab
while read genome 
do 
	cp $output/$genome/phage_number $output/$genome/phage_number_clean
	sed -i "1i $genome" $output/$genome/phage_number_clean
	cat $output/$genome/phage_number_clean | tr '\n' '\t' | cut -f 1,2 > $output/$genome/phage_number_tab
done < $data_path/genome_name

while read genome 
   do 
    accession=`grep "$genome" $data_path/total_genomes/classes.txt | cut -f 2`
    name=`grep "$genome" $data_path/total_genomes/classes.txt | cut -f 3`
    cp $output/$accession.fna/phage_number $output/$accession.fna/phage_number_clean
  	sed -i "1i $name" $output/$accession.fna/phage_number_clean 
 	cat $output/$accession.fna/phage_number_clean | tr '\n' '\t' | cut -f 1,2 > $output/$accession.fna/phage_number_tab
done < $data_path/RefSeq_Accession.tsv

cat $output/*/phage_number_tab | sed s'/.fasta//g' | sed s'/e_Shewanella_putrefaciens_//g' | sed s'/e_//g' | sed s'/_micro12//g' | sed s'/_micro13//g' | sed s'/_micro9//g' | sed s'/_micro22//g' | sed s'/_micro1//g' | sed s'/_1//g' > $output/Total_phage
cat $output/Total_phage | sed s'/ /_/g' | awk '{if ($2>0) {print $0}}' | sed s'/_/ /g' > $output/Total_phage_top