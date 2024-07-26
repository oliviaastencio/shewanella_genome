#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

######################################
#####FUNCTIONS
#####################################
# pendientes (reejecutar:
#tp_clean
# renombrar los resultados de Dfast y de los genomas

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
	grep -v "Assembly Accession" $output/RefSeq.tsv | cut -f 1 > $output/RefSeq_Accession.tsv

	file="$output/RefSeq_Accession.tsv"
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
	
	rm -r $genome_analysis_path/dfast_0000 # $genome_analysis_path/pyani_0000  $genome_analysis_path/Sibelia_0000 $data_path/total_genomes $data_path/all_genome_list
	mkdir -p $data_path/total_genomes
	
	ln -s $data_path/genomes_problem/*fasta $data_path/total_genomes
	ln -s $data_path/genomes_down/genomes/*.fna $data_path/total_genomes

	ls $data_path/genomes_down/genomes > $data_path/genome_candidate_list

	ls $data_path/genomes_problem > $data_path/all_genome_list
	ls $data_path/genomes_down/genomes >> $data_path/all_genome_list
	
	vars=`echo "
		\\$data_path=$project_path'/data',
		\\$ab1_name=$data_path/ab1_name,
		\\$scripts_path=$script_path
		" | tr -d [:space:]`

	AutoFlow -c 1 -s -w $template_path/genome_analysis.af -V $vars -t "1-15:00:00" -o $genome_analysis_path

fi

##################################
## Dfast clean
##################################

if [ "$1" == "dfast_clean" ]; then  # re-name the dfast results file 

	dfast.sh $data_path $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser
	
fi 	

##################################
## TP FLOW
##################################

if [ "$1" == "protein_db" ]; then  

	protein_db.sh  $script_path $data_path $SCRATCH/genome_analysis "Shewanella"
fi 

if [ "$1" == "tp_case" ]; then  

	. ~soft_bio_267/initializes/init_autoflow 
	rm -r $genome_analysis_path/transposon/executions
	 output=$genome_analysis_path/transposon/executions

	while read line 
	do 
		mkdir -p $output/$line
		vars=`echo "
			\\$genome_seq=$data_path/total_genomes/$line,
			\\$prot_database=$data_path/tp_data
			" | tr -d [:space:]`
		AutoFlow -c 1 -s -w $template_path/tpflow -V $vars -o "$output/"$line $2
	done < $data_path/all_genome_list
	
fi 

if [ "$1" == "tp_clean" ]; then  #matrix of disrupted and transposable protein, present or absent. Also re-name Tp_case output files
	Tp.sh $data_path $genome_analysis_path  #sbatch $sh_file_path/
fi

##################################
## TarSynFlow
##################################

if [ "$1" == "genes" ]; then

	output=$genome_analysis_path/genes_identification/Tarsynflow
	mkdir -p $output/proteome

	sbatch $sh_file_path/reduce_prot_redundancy.sh  $data_path $genome_analysis_path

fi

if [ "$1" == "genes_comps" ]; then

	make_all_comps.sh $data_path $script_path $genome_analysis_path $template_path

fi

if [ "$1" == "genes_results_protein" ]; then

	sbatch $sh_file_path/get_all_results.sh $output/results $data_path $output/comps  $script_path get_protein
fi 

if [ "$1" == "genes_results_annotation" ]; then

	get_all_results.sh $output/results $data_path $output/comps  $script_path get_annotations
fi  

if [ "$1" == "seqs" ]; then

	extract_seqs.sh "$output/comps/e_Pdp11_1.fasta/"`head -n 1 data/gen_refs`"/procompart_0000/coord_table_with_strand" $data_path/genomes_problem $output $data_path $script_path 
	
fi 

##################################
## Genomic Islands
##################################


if [ "$1"  == "GI" ]; then 
	HTTP_API_token="1d0a19fa-8c75-4868-7487-92ccadf02b57" 
	input=$genome_analysis_path/dfast_0000/genome_annotation
	output=$genome_analysis_path/genomic_island
	mkdir -p $output
	mkdir -p $output/Island_Viewer4_results

	while read ref_genome
	do 
		genome=`echo $ref_genome | tr ' ' '\t'| cut -f 1` 
		ref_num=`echo $ref_genome | tr ' ' '\t'| cut -f 2` 
		curl -X POST -H '$HTTP_API_token' -Fref_accnum="$ref_num" -Fgenome_file=@$input/$genome'_file'/genome.gbk -Fgenome_name="$genome"_genome -Femail_addr="oliastenciogomez@gmail.com" -Fformat_type="GENBANK" https://www.pathogenomics.sfu.ca/islandviewer/rest/submit/ > $output/Island_Viewer4_results/"$genome"_submit
		grep "token" $output/Island_Viewer4_results/"$genome"_submit | sed s'/    "token": "//g'| sed s'/"//g' > $output/Island_Viewer4_results/"$genome"_jobtoken
	done < $data_path/ref_genome_GI
fi 

if [ "$1" == "GI_result" ]; then 

##Query the status of already submitted genomes. Job_token is a job ID, the number where upload de file. I'ts not the same to my token ID or HTTP API Token
##curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/job_token/ -H 'X-authtoken:your_authentication_token'
	output=$genome_analysis_path/genomic_island/Island_Viewer4_results

	rm -r $output/file_scv
	mkdir -p $output/file_scv
	while read genome 
	do
		jobtoken=`head -n1 $output/"$genome"_jobtoken`
		job_token=`echo $jobtoken`
		curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/$jobtoken/download/csv/ -H '$HTTP_API_token' > $output/"$genome".csv
		mv $output/"$genome".csv $output/file_scv
	done < $data_path/all_genome_list
fi 

if [ "$1" == "GI_clean" ]; then 

	sbatch $sh_file_path/GI.sh $data_path  $genome_analysis_path 
	
fi 

##################################
## Prophage by PHASTEST
##################################

# input the genome file in fasta format 
if [ "$1" == "phage_input" ]; then 

	output=$genome_analysis_path/prophage
	rm $output
	mkdir -p $output

	while read genome 
		do
			mkdir -p $output/$genome
			wget --post-file="$data_path/total_genomes/$genome" "https://phastest.ca/phastest_api" -O $output/$genome/submissions
			cat $output/$genome/submissions | tr ',' '\t' | tr ':' '\t' | sed s'/"//g' | cut -f 2 > $output/$genome/job_id
		done < $data_path/all_genome_list
fi 

# download the PHASTEST results

if [ "$1" == "phage_output" ]; then 

	output=$genome_analysis_path/prophage

	rm $output/*/*.zip $output/*/*json $output/*/*fna $output/*/*txt $output/*/submissions
	while read genome 
	do 
		while read job
		do
			wget "https://phastest.ca/submissions/$job.zip" -O $output/$genome/$job.zip
			cd $output/$genome/
			unzip $job.zip
			cd ../../
		done < $output/$genome/job_id
	grep "region" -c $output/$genome/predicted_phage_regions.json > $output/$genome/phage_number
	done < $data_path/all_genome_list
fi 

if [ "$1" == "phage_clean" ]; then  

	phage.sh $data_path $genome_analysis_path/prophage
	
fi 
## REPORTING
####################################
# /mnt/home/users/pab_001_uma/pedro/html_reporting/template
if [ "$1" == "report" ]; then 
	. ~soft_bio_267/initializes/init_ruby
	rm -r $results_path
	mkdir -p $results_path

	### 16SRNA BLAST results 
	grep -h -v '#' $genome_analysis_path/blastn_0000/blast_16S/blast_* | awk '{if (($3>98.7)&&($4=100)) {print $0}}'| tr ' ' '\t' | sed s'/1_S/\tS/g' | cut -f 1,3,4,16 | sed s'/.fasta//g'| sed s'/_1492R//g' | sed s'/_/ /g'> $results_path/blast_16

	### Complete genome comparison 
	grep -v "Q_Shewanella_Pdp11:126" $genome_analysis_path/pyani_0000/genome_pyani_anim/matrix_identity_1.tab | grep -v "R_scaffold" | cut -f 1,127,128,129,130,131,132,133,134 | sed s'/Shewanella / S./g' | sort >> $results_path/pyani_identity
	sed -i '1ishewanella strains \t Pdp11 \t SH12 \t SH16 \t SH4 \t SH6 \t SH9 \t SdM1 \t SdM2' $results_path/pyani_identity
	grep -v "Q_Shewanella_Pdp11:126" $genome_analysis_path/pyani_0000/genome_pyani_anim/matrix_coverage_1.tab | grep -v "R_scaffold" | cut -f 1,127,128,129,130,131,132,133,134 | sed s'/Shewanella / S./g' | sort >> $results_path/pyani_coverage
	sed -i '1ishewanella strains \t Pdp11 \t SH12 \t SH16 \t SH4 \t SH6 \t SH9 \t SdM1 \t SdM2' $results_path/pyani_coverage 
	
	###### Synteny block and Sibelia 
	cp $genome_analysis_path/Sibelia_0000/e_Pdp11_1/GCF_003052765.1.fna/circos/circos.png $results_path/S_baltica_128:Pdp11.png 
	cp $genome_analysis_path/Sibelia_0000/e_Pdp11_1/GCF_025402875.1.fna/circos/circos.png $results_path/S_putrefaciens_4H:Pdp11.png

	### Genome annotation by DFAST
	cat $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/COG_annotation | sed s'/Shewanella/S./g'> $results_path/COG_annotation
	cat $genome_analysis_path/dfast_0000/genome_annotation/results_dfast_parser/COG_annotation_relative | sed s'/Shewanella/S./g'> $results_path/COG_annotation_relative
	
	### Transposon 
	cat $genome_analysis_path/transposon/executions/Total_tp | sed s'/Shewanella/S./g' > $results_path/Total_tp
	sed -i '1ishewanella strains \t transposable elements' $results_path/Total_tp  
	cat $genome_analysis_path/transposon/executions/Total_tp_relative | sed s'/Shewanella/S./g' > $results_path/Total_tp_relative
	sed -i '1ishewanella strains \t transposable elements' $results_path/Total_tp_relative

	less -S $genome_analysis_path/transposon/executions/e_Pdp11_1.fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,2,3,4 | tr ',' '\t' | cut -f 1,2,3,4 > $genome_analysis_path/transposon/executions/e_Pdp11_1.fasta/transposons_finder.rb_0000/results/Pdp11_tp_1
	less -S $genome_analysis_path/transposon/executions/e_Pdp11_1.fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,5 | tr ',' '\t' | cut -f 1,2 > $genome_analysis_path/transposon/executions/e_Pdp11_1.fasta/transposons_finder.rb_0000/results/Pdp11_tp_2
	merge_tabular.rb $genome_analysis_path/transposon/executions/e_Pdp11_1.fasta/transposons_finder.rb_0000/results/Pdp11_tp_1 $genome_analysis_path/transposon/executions/e_Pdp11_1.fasta/transposons_finder.rb_0000/results/Pdp11_tp_2 > $results_path/Pdp11_tp
	
	### Transposon matrix 
	cp $genome_analysis_path/transposon/executions/all_interrupt_names_tab $results_path/Tp_interrupt
	cp $genome_analysis_path/transposon/executions/all_interrupt_names_tab_total $results_path/Tp_interrupt_total
	cp $genome_analysis_path/transposon/executions/all_transposase_names_tab $results_path/Tp_transposable
	cp $genome_analysis_path/transposon/executions/all_transposase_names_tab_total $results_path/Tp_transposable_total

	### Tarsynflow
	grep -v 'Entry' $genome_analysis_path/genes_identification/Tarsynflow/results/Ref_specific_filtered_annotated_prots | cut -f 1,5,7 > $results_path/specific_genes
	
	### Genomic Island Total
	cat $genome_analysis_path/genomic_island/genomic_island_results/Total_GI | sed s'/Shewanella/S./g' > $results_path/GI_total
	sed -i '1ishewanella strains \t GIs number' $results_path/GI_total

	### Genomic Island COG categories
	grep "categories" $genome_analysis_path/genomic_island/genomic_island_results/Total_category > $genome_analysis_path/genomic_island/genomic_island_results/Total_category_name
	grep -v "categories" $genome_analysis_path/genomic_island/genomic_island_results/Total_category > $genome_analysis_path/genomic_island/genomic_island_results/Total_category_value
	merge_tabular.rb  $data_path/COG_categories_all $genome_analysis_path/genomic_island/genomic_island_results/Total_category_value | cut -f 2-135 > $genome_analysis_path/genomic_island/genomic_island_results/GI_categories
	cat $genome_analysis_path/genomic_island/genomic_island_results/Total_category_name $genome_analysis_path/genomic_island/genomic_island_results/GI_categories | sed s'/Shewanella/S./g' > $results_path/GI_categories

	### Genomic Island COG categories
	grep "categories" $genome_analysis_path/genomic_island/genomic_island_results/Total_category_relative > $genome_analysis_path/genomic_island/genomic_island_results/Total_category_relative_name
	grep -v "categories" $genome_analysis_path/genomic_island/genomic_island_results/Total_category_relative > $genome_analysis_path/genomic_island/genomic_island_results/Total_category_relative_value
	merge_tabular.rb  $data_path/COG_categories_all $genome_analysis_path/genomic_island/genomic_island_results/Total_category_relative_value | cut -f 2-135 > $genome_analysis_path/genomic_island/genomic_island_results/GI_categories_relative
	cat $genome_analysis_path/genomic_island/genomic_island_results/Total_category_relative_name $genome_analysis_path/genomic_island/genomic_island_results/GI_categories_relative | sed s'/Shewanella/S./g' > $results_path/GI_categories_relative

	### Prophage 
	cat $genome_analysis_path/prophage/Total_phage_top| sed s'/Shewanella/S./g' > $results_path/Total_phage
	sed -i '1ishewanella strains \t prophages number' $results_path/Total_phage

	### Report NEW
	paths=`echo -e "
	$results_path/blast_16,
	$results_path/COG_annotation,
	$results_path/COG_annotation_relative,
	$results_path/pyani_identity,
	$results_path/pyani_coverage,
	$results_path/Total_tp,
	$results_path/Total_tp_relative,
	$results_path/Pdp11_tp,
	$results_path/Tp_interrupt,
	$results_path/Tp_interrupt_total,
	$results_path/Tp_transposable,
	$results_path/Tp_transposable_total,
	$results_path/specific_genes,
	$results_path/GI_total,
	$results_path/GI_categories,
	$results_path/GI_categories_relative,
	$results_path/Total_phage
	" | tr -d [:space:]`
	report_html -t $template_path/report_template.erb -d $paths -o $results_path/project_report
	
fi
