#! /usr/bin/env bash

##################################################
###########ab1_CLEAN
##################################################
data_path=$1
source ~soft_bio_267/initializes/init_emboss
module load bbmap/38.50b
module load fastqc

mkdir -p $data_path/16S_file
	main_folder=`pwd`
	while read -r folder_16s
	do
	    mkdir -p $data_path/16S_file/$folder_16s
	    seqret -sformat abi -osformat fastq -auto -stdout -sequence $data_path/ab1'/'$folder_16s'.ab1' > $data_path/16S_file/$folder_16s'/'$folder_16s'.fastq'
	    bbduk.sh -Xmx1g in=$data_path/16S_file/$folder_16s'/'$folder_16s'.fastq' out=$data_path/16S_file/$folder_16s'/'$folder_16s'clean.fastq' qtrim=rl trimq=20 qin=33
	    seqret -sformat fastq -osformat fasta -auto -stdout -sequence $data_path/16S_file/$folder_16s'/'$folder_16s'clean.fastq' > $data_path/16S_file/$folder_16s'/'$folder_16s'.fasta'
	    mkdir -p $data_path/16S_file/$folder_16s'/'analysis
	    fastqc -o $data_path/16S_file/$folder_16s'/'analysis $data_path/16S_file/$folder_16s'/'$folder_16s'.fastq'   
	done < $data_path/ab1_name

######################################################
######################OPTIONAL STEP###################
######################################################

sed s'/JM_1492R/Pdp11_1492R/g' $data_path/16S_file/Pdp11_1492R/Pdp11_1492R.fasta > $data_path/16S_file/Pdp11_1492R/Pdp11_1492R_1.fasta
rm $data_path/16S_file/Pdp11_1492R/Pdp11_1492R.fasta
sed s'/FMCCB65_1492R/SdM1_1492R/g' $data_path/16S_file/SdM1_1492R/SdM1_1492R.fasta > $data_path/16S_file/SdM1_1492R/SdM1_1492R_1.fasta
rm $data_path/16S_file/SdM1_1492R/SdM1_1492R.fasta
sed s'/CCCT_1492R/SdM2_1492R/g' $data_path/16S_file/SdM2_1492R/SdM2_1492R.fasta > $data_path/16S_file/SdM2_1492R/SdM2_1492R_1.fasta
rm $data_path/16S_file/SdM2_1492R/SdM2_1492R.fasta
sed s'/SH2_1492R/SH12_1492R/g' $data_path/16S_file/SH12_1492R/SH12_1492R.fasta > $data_path/16S_file/SH12_1492R/SH12_1492R_1.fasta
rm $data_path/16S_file/SH12_1492R/SH12_1492R.fasta
sed s'/SH8_1492R/SH6_1492R.fasta/g' $data_path/16S_file/SH6_1492R/SH6_1492R.fasta > $data_path/16S_file/SH6_1492R/SH6_1492R_1.fasta
rm $data_path/16S_file/SH6_1492R/SH6_1492R.fasta
sed s'/SH19_1492R/SH9_1492R/g' $data_path/16S_file/SH9_1492R/SH9_1492R.fasta > $data_path/16S_file/SH9_1492R/SH9_1492R_1.fasta
rm $data_path/16S_file/SH9_1492R/SH9_1492R.fasta
cp $data_path/16S_file/SH16_1492R/SH16_1492R.fasta $data_path/16S_file/SH16_1492R/SH16_1492R_1.fasta
cp $data_path/16S_file/SH4_1492R/SH4_1492R.fasta $data_path/16S_file/SH4_1492R/SH4_1492R_1.fasta