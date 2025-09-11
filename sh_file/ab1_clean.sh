#! /usr/bin/env bash

##################################################
###########ab1_CLEAN
##################################################
data_path=$1

source ~soft_bio_267/initializes/init_emboss
module load bbmap/38.50b
module load fastqc/0.11.9

rm -r $data_path/16S_file
mkdir -p $data_path/16S_file
	main_folder=`pwd`
	while read -r folder_16s
	do
	    mkdir -p $data_path/16S_file/$folder_16s
	    seqret -sformat abi -osformat fastq -auto -stdout -sequence $data_path/ab1'/'$folder_16s'.ab1' > $data_path/16S_file/$folder_16s'/'$folder_16s'.fastq'
	    bbduk.sh -Xmx1g in=$data_path/16S_file/$folder_16s'/'$folder_16s'.fastq' out=$data_path/16S_file/$folder_16s'/'$folder_16s'clean.fastq' qtrim=rl trimq=20 qin=33
	    seqret -sformat fastq -osformat fasta -auto -stdout -sequence $data_path/16S_file/$folder_16s'/'$folder_16s'clean.fastq' | grep -v ">"  > $data_path/16S_file/$folder_16s'/'$folder_16s'.fasta' 
	    echo -e ">$folder_16s" | cat - $data_path/16S_file/$folder_16s'/'$folder_16s'.fasta' > temp && mv temp $data_path/16S_file/$folder_16s'/'$folder_16s'.fasta'
	    mkdir -p $data_path/16S_file/$folder_16s'/'analysis
	    fastqc -o $data_path/16S_file/$folder_16s'/'analysis $data_path/16S_file/$folder_16s'/'$folder_16s'.fastq'   
	done < $data_path/ab1_name
