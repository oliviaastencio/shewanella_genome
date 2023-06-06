#! /usr/bin/env bash

PATH=~soft_bio_267/programs/x86_64/scripts/:$PATH
summary=$1  
protein_fasta=$2  
out_put=$3  

rm -r $out_put/comparative
mkdir -p $out_put/comparative
module load blast_plus/2.2.30+
module load cdhit/4.5.4
awk '{ print $2":"$3 }' $summary > $out_put/comparative/tranposase_coordinates_length

################################################
#################transposase sequence extraction
################################################
while read sequence
do 
	folder='tp_'`echo $sequence | tr ':' '_'`
	mkdir -p $out_put/comparative/$folder

	while read coordinate
	do 
		echo '>'$folder'_'$coordinate > $out_put/comparative/$folder/transposase_"$coordinate".fasta
		length=`echo $coordinate | tr ':' '-'`
		head -n2 $out_put/e_Pdp11_1/lista_to_fasta.rb_0000/tp_case/$folder/tp.fasta | cut -c $length | awk 'NF' >> $out_put/comparative/$folder/transposase_"$coordinate".fasta
	done < $out_put/comparative/tranposase_coordinates_length
	cat $out_put/comparative/$folder/transposase_* > $out_put/comparative/$folder/transposase.fasta
	grep $sequence $summary | cut -f 4 | tr ',' '\n' | tr '\t' '\n' > $out_put/comparative/$folder/protein_code
	lista_to_fasta.rb $protein_fasta $out_put/comparative/$folder/protein_code > $out_put/comparative/$folder/protein.fasta 
	cd-hit -i $out_put/comparative/$folder/protein.fasta -o $out_put/comparative/$folder/protein_new.fasta -c 0.95	
done <$out_put/e_Pdp11_1/lista_to_fasta.rb_0000/tp_case/tp_transposons_finder_coordinates

#################################################
############transposase comparative##########
###############################################

mkdir -p $out_put/comparative/blastx_comparative
cat $out_put/comparative/tp_*/transposase.fasta > $out_put/comparative/total_Pdp11_transposase.fasta
blastn -query $out_put/comparative/total_Pdp11_transposase.fasta -subject $out_put/comparative/total_Pdp11_transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_total
grep -v "#" $out_put/comparative/blastx_comparative/blast_total | awk '{if ($3>=35) {print $0}}' | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary
