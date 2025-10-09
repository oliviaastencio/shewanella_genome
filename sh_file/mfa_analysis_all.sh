#!/usr/bin/env bash
##JOB_GROUP_ID=template_CA.txt_1736508105
#SBATCH --cpus-per-task=1
#SBATCH --mem='80gb'
#SBATCH --time='7-00:00:00'
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

data=$1
output=$2

source ~soft_bio_267/initializes/init_degenes_hunter
source ~soft_bio_267/initializes/init_R  
mkdir results
multivar_mine.R -i "$data/GI_Functions,$data/Whole_Genome_Functions,$data/Relative_Mobilome_Number" -A "GI_Functions:s,Whole_Genome_Functions:s,Relative_Mobilome_Number:s" -c -o results/$output --seed 123
