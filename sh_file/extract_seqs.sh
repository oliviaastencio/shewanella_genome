#! /usr/bin/env bash

coord_table=$1 #"genes_identification/Tarsynflow/comps/e_Pdp11_1.fasta/"`head -n 1 data/gen_refs_1`"/procompart_0000/coord_table_with_strand"
genome_path=$2
output=$3
data=$4
script=$5
genome="$2/e_Shewanella_Pdp11.fasta"
coord_table_path=$3/comps/Shewanella_Pdp11.fasta #"genes_identification/Tarsynflow/comps/e_Pdp11_1.fasta"


mkdir -p $output/extracted_seqs
cut -f 1,2,5 $output/results/Query_specific_filtered_annotated_prots | grep -v -i fragment > $output/extracted_seqs/seqs_list_query.txt
cut -f 1  $output/extracted_seqs/seqs_list_query.txt > $output/extracted_seqs/seqs_query.lst
grep -F -f $output/extracted_seqs/seqs_query.lst $coord_table > $output/extracted_seqs/seqs_coords_query.txt
$script/scripts_Tarsynflow/get_seqs.rb -c $output/extracted_seqs/seqs_coords_query.txt -f $genome -r 150 > $output/extracted_seqs/seqs_query.fasta


cut -f 1,2,5 $output/results/Ref_specific_filtered_annotated_prots | grep -v -i fragment > $output/extracted_seqs/seqs_list_ref.txt
cut -f 1  $output/extracted_seqs/seqs_list_ref.txt > $output/extracted_seqs/seqs_ref.lst

while read REF
do
        grep -F -f $output/extracted_seqs/seqs_ref.lst $coord_table_path"/"$REF"/procompart_0001/coord_table_with_strand" > $output/extracted_seqs/seqs_coords_ref_"$REF".txt
        $script/scripts_Tarsynflow/get_seqs.rb -c $output/extracted_seqs/seqs_coords_ref_"$REF".txt -f $genome_path"/e_"$REF -r 150 > $output/extracted_seqs/seqs_ref_"$REF".fasta
done < $data/gen_refs

