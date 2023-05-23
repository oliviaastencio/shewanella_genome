#! /usr/bin/env bash

project_path=`pwd`
data_path=$project_path'/data'
mkdir -p $data_path/16S_Shewanella
	main_folder=`pwd`
	while read -r line
	do
	    folder_16s=`echo $line`
	    mkdir -p $data_path/16S_Shewanella/$folder_16s
	    source ~soft_bio_267/initializes/init_emboss
	    seqret -sformat abi -osformat fastq -auto -stdout -sequence $data_path/ab1'/'$folder_16s'.ab1' > $data_path/16S_Shewanella/$folder_16s'/'$folder_16s'.fastq'
	    module load bbmap/38.50b
	    bbduk.sh -Xmx1g in=$data_path/16S_Shewanella/$folder_16s'/'$folder_16s'.fastq' out=$data_path/16S_Shewanella/$folder_16s'/'$folder_16s'clean.fastq' qtrim=rl trimq=20 qin=33
	    seqret -sformat fastq -osformat fasta -auto -stdout -sequence $data_path/16S_Shewanella/$folder_16s'/'$folder_16s'clean.fastq' > $data_path/16S_Shewanella/$folder_16s'/'$folder_16s'.fasta'
	    mkdir -p $data_path/16S_Shewanella/$folder_16s'/'analysis
	    module load fastqc
	    fastqc -o $data_path/16S_Shewanella/$folder_16s'/'analysis $data_path/16S_Shewanella/$folder_16s'/'$folder_16s'.fastq'   
	done < $data_path/ab1_name

	sed s'/JM_1492R/Pdp11_1492R/g' $data_path/16S_Shewanella/Pdp11_1492R/Pdp11_1492R.fasta > $data_path/16S_Shewanella/Pdp11_1492R/Pdp11_1492R_1.fasta
	rm $data_path/16S_Shewanella/Pdp11_1492R/Pdp11_1492R.fasta
	sed s'/FMCCB65_1492R/SdM1_1492R/g' $data_path/16S_Shewanella/SdM1_1492R/SdM1_1492R.fasta > $data_path/16S_Shewanella/SdM1_1492R/SdM1_1492R_1.fasta
	rm $data_path/16S_Shewanella/SdM1_1492R/SdM1_1492R.fasta
	sed s'/CCCT_1492R/SdM2_1492R/g' $data_path/16S_Shewanella/SdM2_1492R/SdM2_1492R.fasta > $data_path/16S_Shewanella/SdM2_1492R/SdM2_1492R_1.fasta
	rm $data_path/16S_Shewanella/SdM2_1492R/SdM2_1492R.fasta
	sed s'/SH2_1492R/SH12_1492R/g' $data_path/16S_Shewanella/SH12_1492R/SH12_1492R.fasta > $data_path/16S_Shewanella/SH12_1492R/SH12_1492R_1.fasta
	rm $data_path/16S_Shewanella/SH12_1492R/SH12_1492R.fasta
	sed s'/SH8_1492R/SH6_1492R.fasta/g' $data_path/16S_Shewanella/SH6_1492R/SH6_1492R.fasta > $data_path/16S_Shewanella/SH6_1492R/SH6_1492R_1.fasta
	rm $data_path/16S_Shewanella/SH6_1492R/SH6_1492R.fasta
	sed s'/SH19_1492R/SH9_1492R/g' $data_path/16S_Shewanella/SH9_1492R/SH9_1492R.fasta > $data_path/16S_Shewanella/SH9_1492R/SH9_1492R_1.fasta
	rm $data_path/16S_Shewanella/SH9_1492R/SH9_1492R.fasta
	cp $data_path/16S_Shewanella/SH16_1492R/SH16_1492R.fasta $data_path/16S_Shewanella/SH16_1492R/SH16_1492R_1.fasta
	cp $data_path/16S_Shewanella/SH4_1492R/SH4_1492R.fasta $data_path/16S_Shewanella/SH4_1492R/SH4_1492R_1.fasta