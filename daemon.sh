#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

######################################
#####FUNCTIONS
#####################################

##################################################
######Complete genomes database
#################################################
function genomes_database {
	output=$1
	module load ruby/2.4.1 
	rm -r $output/genomes_down
	mkdir -p $output/genomes_down/sequence_download

	input=$output/genomes_down/sequence_download 
	
	# #######We online need the genomic.fna.gz 
	grep -v "Assembly Accession" $output/RefSeq.tsv | cut -f 1 > $output/RefSeq_list.tsv

	file="$output/RefSeq_list.tsv"
	while IFS= read -r line
	do
	    wget "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$line/download?include_annotation_type=GENOME_FASTA" -O genome.zip
	    unzip genome.zip
	  	mv ncbi_dataset/data/$line/*.fna $output/genomes_down/sequence_download/$line.fna
	    rm -r genome.zip ncbi_dataset README.md
	done < "$file" 

	sed s'/plasmid/plasmid complete genome/g'  $output/genomes_down/sequence_download/GCF_949794895.1.fna | sed s'/chromosome/complete/g' > $output/genomes_down/sequence_download/GCF_949794895.1_2.fna
	sed s'/plasmid/plasmid complete genome/g'  $output/genomes_down/sequence_download/GCF_963676895.1.fna | sed s'/chromosome/complete/g' > $output/genomes_down/sequence_download/GCF_963676895.1_2.fna
	sed s'/plasmid/plasmid complete genome/g'  $output/genomes_down/sequence_download/GCF_024363255.1.fna | sed s'/chromosome/complete/g' > $output/genomes_down/sequence_download/GCF_024363255.1_2.fna
	rm $output/genomes_down/sequence_download/GCF_949794895.1.fna $output/genomes_down/sequence_download/GCF_024363255.1.fna $output/genomes_down/sequence_download/GCF_963676895.1.fna
	######WE Separated genomes and plasmids
	genome_folder=$output/genomes_down/genomes
	plasmid_folder=$output/genomes_down/plasmids
	mkdir -p $genome_folder $plasmid_folder

	strains=( `ls $input` )
	for strain in "${strains[@]}"
	do
		genome_parser.rb $input/$strain $genome_folder/$strain $plasmid_folder/$strain
	done
}

#####################################################
#######  MAIN 
##################################################
########1- genome reassigment (16S RNA, complete genome comparison and genome anotation) ####################  

project_path=`pwd`
data_path=$project_path'/data'
script_path=$project_path'/script'
sh_file_path=$project_path'/sh_file'
template_path=$project_path'/templates'
results_path=$project_path'/results'
export PATH=$script_path:$PATH
export PATH=$script_path'/scripts_Tarsynflow:'$PATH
export PATH=$sh_file_path:$PATH

genome_analysis_path=$SCRATCH/genome_analysis
transposon_analysis_path=$SCRATCH/transposon_result

if [ "$1" == "down" ]; then  

	genomes_database $data_path
fi
 	
if [ "$1" == "ab1_clean" ]; then   
	ab1_clean.sh $data_path
fi

if [ "$1" == "char" ]; then  
	. ~soft_bio_267/initializes/init_autoflow  

	rm -r $genome_analysis_path/pyani_0000  $genome_analysis_path/Sibelia_0000 $data_path/total_genomes $data_path/all_genome_list
	mkdir -p $data_path/total_genomes
	
	ln -s $data_path/genomes_problem/*fasta $data_path/total_genomes
	ln -s $data_path/genomes_down/genomes/*.fna $data_path/total_genomes

	ls $data_path/genomes_down/genomes > $data_path/genome_candidate_list

	ls $data_path/genomes_problem > $data_path/all_genome_list
	ls $data_path/genomes_down/genomes >> $data_path/all_genome_list
	
	vars=`echo "
		\\$data_path=$project_path'/data',
		\\$ab1_name=$data_path/ab1_name,
		\\$genome_list=$data_path/genomes_link,
		\\$scripts_path=$script_path
		" | tr -d [:space:]`

	AutoFlow -c 1 -s -w $template_path/genome_analysis.af -V $vars -t "1-10:00:00" -o $genome_analysis_path
fi


##################################
## TP FLOW
##################################

if [ "$1" == "protein_db" ]; then  

	protein_db.sh  $script_path $data_path $SCRATCH/genome_analysis "Shewanella"
fi 

if [ "$1" == "tp_case" ]; then  
	
	. ~soft_bio_267/initializes/init_autoflow 
	while read line 
	do 
		rm -r $genome_analysis_path/transposon/executions/$line
		mkdir -p $genome_analysis_path/transposon/executions/$line
		vars=`echo "
			\\$genome_seq=$data_path/genomes_problem/$line.fasta,
			\\$prot_database=$data_path/tp_data
			" | tr -d [:space:]`
		AutoFlow -c 1 -s -w $template_path/tpflow -V $vars -o "$genome_analysis_path/transposon/executions/"$line $2
	done < $data_path/genome_name
fi 

if [ "$1" == "tp_comparative" ]; then   #####OPTIONAL STEP##### LISTO 

	tp_comparative.sh $genome_analysis_path/transposon/executions/e_Pdp11_1/transposons_finder.rb_0000/results/summary.txt $data_path/tp_data/total_prots.fasta $genome_analysis_path/transposon/executions

fi 

##################################
## TarSynFlow
##################################

if [ "$1" == "genes" ]; then

	mkdir -p $genome_analysis_path/genes_identification/Tarsynflow/proteome

	sbatch $sh_file_path/reduce_prot_redundancy.sh  $data_path $genome_analysis_path

fi

if [ "$1" == "genes_comps" ]; then

	make_all_comps.sh $data_path $script_path $genome_analysis_path $template_path

fi

if [ "$1" == "genes_results" ]; then

	sbatch $sh_file_path/get_all_results.sh $genome_analysis_path/genes_identification/Tarsynflow/results $data_path $genome_analysis_path/genes_identification/Tarsynflow/comps  $script_path 

fi 

if [ "$1" == "seqs" ]; then

	extract_seqs.sh "$genome_analysis_path/genes_identification/Tarsynflow/comps/e_Pdp11_1.fasta/"`head -n 1 data/gen_refs`"/procompart_0000/coord_table_with_strand" $data_path/genomes_problem $genome_analysis_path/genes_identification/Tarsynflow $data_path $script_path 
	
fi 

##################################
## Genomic Islands
##################################
HTTP_API_token="1d0a19fa-8c75-4868-7487-92ccadf02b57" 
input=$genome_analysis_path/dfast_0000/genome_annotation
output=$genome_analysis_path/genomic_island
mkdir -p $output
mkdir -p $output/Island_Viewer4_results

if [ "$1"  == "GI" ]; then 
	while read ref_genome
	do 
		genome=`echo $ref_genome | tr ' ' '\t'| cut -f 1` 
		ref_num=`echo $ref_genome | tr ' ' '\t'| cut -f 2` 
		curl -X POST -H '$HTTP_API_token' -Fref_accnum="$ref_num" -Fgenome_file=@$input/$genome/genome.gbk -Fgenome_name="$genome"_genome -Femail_addr="oliastenciogomez@gmail.com" -Fformat_type="GENBANK" https://www.pathogenomics.sfu.ca/islandviewer/rest/submit/ > $output/Island_Viewer4_results/"$genome"_submit
		grep "token" $output/Island_Viewer4_results/"$genome"_submit | sed s'/    "token": "//g'| sed s'/"//g' > $output/Island_Viewer4_results/"$genome"_jobtoken
	done < $data_path/ref_genome_GI
fi 

if [ "$1" == "GI_result" ]; then 

##Query the status of already submitted genomes. Job_token is a job ID, the number where upload de file. I'ts not the same to my token ID or HTTP API Token
##curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/job_token/ -H 'X-authtoken:your_authentication_token'

	mkdir -p $genome_analysis_path/genomic_island/Island_Viewer4_results/file_scv
	while read genome 
	do
		jobtoken=`head -n1 $genome_analysis_path/genomic_island/Island_Viewer4_results/"$genome"_jobtoken`
		job_token=`echo $jobtoken`
		curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/$jobtoken/download/csv/ -H '$HTTP_API_token' > $genome_analysis_path/genomic_island/Island_Viewer4_results/"$genome".csv
		mv $genome_analysis_path/genomic_island/Island_Viewer4_results/"$genome".csv $genome_analysis_path/genomic_island/Island_Viewer4_results/file_scv
	done < $data_path/genome_name
fi 

if [ "$1" == "GI_filtre" ]; then 

	module load ruby/2.4.1
	input=$genome_analysis_path/genomic_island/Island_Viewer4_results/file_scv 
	output=$genome_analysis_path/genomic_island/genomic_island_results
	
	rm -r $output
	mkdir -p $output
	strains_island=( `ls $input` )
	for strain in "${strains_island[@]}"
	do
		mkdir -p $output/$strain
		island_filtre.rb $input/$strain $output/$strain/$strain'_length'
		grep 'Predicted by at least' $output/$strain/$strain'_length'|cut -f 1,2,3,6,7,11,12 > $output/$strain/$strain'_Integrated'
		sed -i "1i Island_start\tIsland_end\tLength\tGene_ID\tLocus\tProduct\tExternal_Annotations" $output/$strain/$strain'_Integrated'
		echo $strain > $output/$strain/Total_genomic_island
		grep -v 'Island_start' $output/$strain/$strain'_Integrated' | cut -f1 | sort -u | wc -l >> $output/$strain/Total_genomic_island
		grep -v " " $output/$strain/Total_genomic_island | tr '\n' '\t' >> $output/$strain/Total_genomic_island_tab
	done 
	cat $output/*/*tab | sed s'/e_S/\nS/g' | sed s'/.csv//g' > $output/Total_GI
	
fi 

## REPORTING
####################################
# /mnt/home/users/pab_001_uma/pedro/html_reporting/template
if [ "$1" == "report" ]; then 
	. ~soft_bio_267/initializes/init_ruby
	rm -r $results_path
	mkdir -p $results_path

	# Get tabular data
	grep -h -v '#' $genome_analysis_path/blastn_0000/blast_16S/blast_* | cut -f 1,2,3,15 > $results_path/blast_16
	cp $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/Total_table.txt $results_path/COG_annotation
	merge_tabular.rb  $data_path/COG_categories $results_path/COG_annotation | grep -v "category" | cut -f 2,3,4,5,6,7,8,9,10 > $results_path/COG_annotation_complete
	sed -i '1icategory \t Pdp11 \t SdM1 \t SdM2 \t SH12 \t SH16 \t SH4 \t SH6 \t SH9' $results_path/COG_annotation_complete
	
	less -S $genome_analysis_path/transposon/executions/e_Pdp11_1/transposons_finder.rb_0000/results/summary.txt | cut -f 1,2,3,4 | tr ',' '\t' | cut -f 1,2,3,4 > $results_path/Pdp11_tp_1
	less -S $genome_analysis_path/transposon/executions/e_Pdp11_1/transposons_finder.rb_0000/results/summary.txt | cut -f 1,5 | tr ',' '\t' | cut -f 1,2 > $results_path/Pdp11_tp_2
	merge_tabular.rb $results_path/Pdp11_tp_1 $results_path/Pdp11_tp_2 > $results_path/Pdp11_tp

	cat $genome_analysis_path/transposon/executions/comparative/blastx_comparative/blast_summary | tr '_' '\t' | cut -f 2,6,9,10,16 | sort -u > $results_path/tp_comparative
	grep -v 'Entry' $genome_analysis_path/genes_identification/Tarsynflow/results/Ref_specific_filtered_annotated_prots | cut -f 1,5,7 > $results_path/specific_genes
	
	####### verify de id number for our genomes problems (120-127), that correspond to column (121-128)
	
	grep "Shewanella " $genome_analysis_path/pyani_0000/genome_pyani_anim/matrix_identity_1.tab | cut -f 1,121,122,123,124,125,126,127,128 | sed s'/Shewanella /S./g' | sort >> $results_path/pyani_identity
	sed -i '1ishewanella strains \t Pdp11 \t SH12 \t SH16 \t SH4 \t SH6 \t SH9 \t SdM1 \t SdM2' $results_path/pyani_identity
	grep "Shewanella " $genome_analysis_path/pyani_0000/genome_pyani_anim/matrix_coverage_1.tab | cut -f 1,121,122,123,124,125,126,127,128 | sed s'/Shewanella /S./g' | sort >> $results_path/pyani_coverage
	sed -i '1ishewanella strains \t Pdp11 \t SH12 \t SH16 \t SH4 \t SH6 \t SH9 \t SdM1 \t SdM2' $results_path/pyani_coverage
	###### $data_path/total_genomes/classes.txt was used to verify the assembly number to strain name

	cp $genome_analysis_path/Sibelia_0000/e_Pdp11_1/GCF_003052765.1_ASM305276v1_genomic.fna/circos/circos.png $results_path/S_baltica_128:Pdp11.png 
	cp $genome_analysis_path/Sibelia_0000/e_Pdp11_1/GCF_025402875.1_ASM2540287v1_genomic.fna/circos/circos.png $results_path/S_putrefaciens_4H:Pdp11.png
	
	cp $genome_analysis_path/genomic_island/genomic_island_results/Total_GI $results_path/GI_total
	sed -i '1ishewanella strains \t GIs number' $results_path/GI_total
	cut -f 1,2,3,4,5,6 $genome_analysis_path/genomic_island/genomic_island_results/e_Pdp11_1.csv/e_Pdp11_1.csv_Integrated > $results_path/GI_Pdp11

	paths=`echo -e "
	$results_path/blast_16,
	$results_path/COG_annotation_complete,
	$results_path/pyani_identity,
	$results_path/pyani_coverage,
	$results_path/Pdp11_tp,
	$results_path/tp_comparative,
	$results_path/specific_genes,
	$results_path/GI_total,
	$results_path/GI_Pdp11
	" | tr -d [:space:]`
	report_html -t $template_path/report_template.erb -d $paths -o $results_path/project_report

	###Nota:pendiente subir imagenes del Sibelia 

fi
