#! /usr/bin/env bash

PATH=~soft_bio_267/programs/x86_64/scripts/:$PATH
summary=/mnt/scratch/users/pab_001_uma/oliastencio/genome_analysis/transposon/executions/e_Pdp11_1/transposons_finder.rb_0000/results/summary.txt  #shewanella_genomes_paper/data_script/transposons/tp_finder_results/results/summary.txt
protein_fasta=data/tp_data/total_prots.fasta   #shewanella_genomes_paper/data_script/transposons/data/total_prots.fasta
out_put=/mnt/scratch/users/pab_001_uma/oliastencio/genome_analysis/transposon/executions   #shewanella_genomes_paper/result/transposon/comparative

	#transposon complete sequence = disrupted protein fragment + transposon sequence + disrupted protein fragment
mkdir -p $out_put/comparative
module load blast_plus/2.2.30+
module load cdhit/4.5.4
awk '{ print $2":"$3 }' $summary > $out_put/comparative/tranposase_coordinates_length

	main_folder=`pwd`

	while read sequence
	do 
		folder='tp_'`echo $sequence | tr ':' '_'`
		call_folder=`echo $folder`
		mkdir -p $out_put/comparative/$folder

			while read coordinate
			do 
				echo '>'$folder'_'$coordinate > $out_put/comparative/$folder/transposase_"$coordinate".fasta
				length=`echo $coordinate | tr ':' '-'`
				head -n2 $out_put/e_Pdp11_1/lista_to_fasta.rb_0000/tp_case/$call_folder/tp.fasta | cut -c $length | awk 'NF' >> $out_put/comparative/$folder/transposase_"$coordinate".fasta
			done < $out_put/comparative/tranposase_coordinates_length
			cat $out_put/comparative/$folder/transposase_* > $out_put/comparative/$folder/transposase.fasta
			grep $sequence $summary | cut -f 4 | tr ',' '\n' | tr '\t' '\n' > $out_put/comparative/$folder/protein_code
			lista_to_fasta.rb $protein_fasta $out_put/comparative/$folder/protein_code > $out_put/comparative/$folder/protein.fasta
			#We have used the transposon proteins sequence. 
			cd-hit -i $out_put/comparative/$folder/protein.fasta -o $out_put/comparative/$folder/protein_new.fasta -c 0.95

	done <$out_put/e_Pdp11_1/lista_to_fasta.rb_0000/tp_case/tp_transposons_finder_coordinates

mkdir -p $out_put/comparative/blastx_comparative
		####comparativa de las 3 parejas de tranposones similares por coordenadas.... tp_850951_856278 y tp_850951_856278; tp_1906106_1911827 y tp_1906106_1911827; también incluimos una tercera pareja aunque no tienen resultados similares 
blastn -query $out_put/comparative/tp_850951_856278/transposase.fasta -subject $out_put/comparative/tp_1261471_1266798/transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_1261471_1266798
awk '{if ($3>=35) {print $0}}' $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_1261471_1266798 | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary_tp_850951_856278_tp_1261471_1266798
blastn -query $out_put/comparative/tp_1906106_1911827/transposase.fasta -subject $out_put/comparative/tp_3073685_3079406/transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_tp_1906106_1911827_tp_3073685_3079406
awk '{if ($3>=35) {print $0}}' $out_put/comparative/blastx_comparative/blast_tp_1906106_1911827_tp_3073685_3079406 | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary_tp_1906106_1911827_tp_3073685_3079406
blastn -query $out_put/comparative/tp_2693976_2700579/transposase.fasta -subject $out_put/comparative/tp_4349020_4355949/transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_tp_2693976_2700579_tp_4349020_4355949
awk '{if ($3>=35) {print $0}}' $out_put/comparative/blastx_comparative/blast_tp_2693976_2700579_tp_4349020_4355949 | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary_tp_2693976_2700579_tp_4349020_4355949

		#####comparativa entre diferentes grupos
blastn -query $out_put/comparative/tp_850951_856278/transposase.fasta -subject $out_put/comparative/tp_1906106_1911827/transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_1906106_1911827
awk '{if ($3>=35) {print $0}}' $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_1906106_1911827 | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary_tp_850951_856278_tp_1906106_1911827
blastn -query $out_put/comparative/tp_850951_856278/transposase.fasta  -subject $out_put/comparative/tp_2693976_2700579/transposase.fasta -outfmt '7 std slen qlen' > $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_2693976_2700579
awk '{if ($3>=35) {print $0}}' $out_put/comparative/blastx_comparative/blast_tp_850951_856278_tp_2693976_2700579 | cut -f 1,2,3,4,6,7,8,13,14 > $out_put/comparative/blastx_comparative/blast_summary_tp_850951_856278_tp_2693976_2700579
