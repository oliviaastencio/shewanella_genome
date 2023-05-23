#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=20:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --constraint=cal

project_path=`pwd`
data_path=$project_path'/data'
genome_analysis_path=$SCRATCH/genome_analysis

hostname
module load cdhit
before=`grep -c '>' $data_path/tp_data/total_prots.fasta`
time cd-hit -i $data_path/tp_data/total_prots.fasta -o $genome_analysis_path/genes_identification/Tarsynflow/proteome/prots_clean.fasta -c 0.6 -M 0 -n 4
after=`grep -c '>' $genome_analysis_path/genes_identification/Tarsynflow/proteome/prots_clean.fasta`
echo -e "Downloaded protein\t$before"
echo -e "Final count\t$after"
echo -e "Percentage\t`echo -e "scale=2; $after*100/$before" | bc`"
