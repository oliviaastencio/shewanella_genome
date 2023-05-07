#! /usr/bin/env bash

coord_table="genes_identification/Tarsynflow/comps/e_Pdp11_1.fasta/"`head -n 1 data/gen_refs_1`"/procompart_0000/coord_table_with_strand"
genome="data/genomes_problem/e_e_Pdp11_1.fasta"
genome_path="data/genomes_problem"
coord_table_path="genes_identification/Tarsynflow/comps/e_Pdp11_1.fasta"

mkdir -p genes_identification/Tarsynflow/extracted_seqs
cut -f 1,2,5 genes_identification/Tarsynflow/results/Query_specific_filtered_annotated_prots | grep -v -i fragment > genes_identification/Tarsynflow/extracted_seqs/seqs_list_query.txt
cut -f 1  genes_identification/Tarsynflow/extracted_seqs/seqs_list_query.txt > genes_identification/Tarsynflow/extracted_seqs/seqs_query.lst
grep -F -f genes_identification/Tarsynflow/extracted_seqs/seqs_query.lst $coord_table > genes_identification/Tarsynflow/extracted_seqs/seqs_coords_query.txt
script/scripts_Tarsynflow/get_seqs.rb -c genes_identification/Tarsynflow/extracted_seqs/seqs_coords_query.txt -f $genome -r 150 > genes_identification/Tarsynflow/extracted_seqs/seqs_query.fasta


cut -f 1,2,5 genes_identification/Tarsynflow/results/Ref_specific_filtered_annotated_prots | grep -v -i fragment > genes_identification/Tarsynflow/extracted_seqs/seqs_list_ref.txt
cut -f 1  genes_identification/Tarsynflow/extracted_seqs/seqs_list_ref.txt > genes_identification/Tarsynflow/extracted_seqs/seqs_ref.lst

while read REF
do
        grep -F -f genes_identification/Tarsynflow/extracted_seqs/seqs_ref.lst $coord_table_path"/"$REF"/procompart_0001/coord_table_with_strand" > genes_identification/Tarsynflow/extracted_seqs/seqs_coords_ref_"$REF".txt
        script/scripts_Tarsynflow/get_seqs.rb -c genes_identification/Tarsynflow/extracted_seqs/seqs_coords_ref_"$REF".txt -f $genome_path"/e_"$REF -r 150 > genes_identification/Tarsynflow/extracted_seqs/seqs_ref_"$REF".fasta
done < data/gen_refs_1

