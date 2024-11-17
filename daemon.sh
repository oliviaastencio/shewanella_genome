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
	
	######## genome download from 

	while read line
	do
		assembly=`grep "$line" $output/RefSeq.tsv | cut -f 1`
		name=`grep "$line" $output/RefSeq.tsv | cut -f 3,5 | tr '\t' '_' | tr ' ' '_'`   # | tr '.' '_' | tr '-' '_'
	    wget "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$assembly/download?include_annotation_type=GENOME_FASTA" -O genome.zip
	    unzip genome.zip
	  	mv ncbi_dataset/data/$assembly/*.fna $output/genomes_down/sequence_download/$name.fasta
	    rm -r genome.zip ncbi_dataset README.md  md5sum.txt
	done < $output/RefSeq.tsv

	sed s'/plasmid/plasmid complete genome/g'  $output/genomes_down/sequence_download/Shewanella_baltica_SF1039.fasta | sed s'/chromosome/complete/g' > $output/genomes_down/sequence_download/Shewanella_baltica_SF1039_1.fasta
	sed s'/plasmid/plasmid complete genome/g'  $output/genomes_down/sequence_download/uncultured_Shewanella_sp._.fasta | sed s'/chromosome/complete/g' > $output/genomes_down/sequence_download/uncultured_Shewanella_sp.fasta
	sed s'/plasmid/plasmid complete genome/g'  $output/genomes_down/sequence_download/Shewanella_sp._NFH-SH190041_.fasta | sed s'/chromosome/complete/g' > $output/genomes_down/sequence_download/Shewanella_sp._NFH-SH190041.fasta
	rm $output/genomes_down/sequence_download/Shewanella_baltica_SF1039.fasta $output/genomes_down/sequence_download/uncultured_Shewanella_sp._.fasta $output/genomes_down/sequence_download/Shewanella_sp._NFH-SH190041_.fasta
	
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

if [ "$1" == "down" ]; then  ######LISTO

	genomes_database $data_path   ######LISTO
fi
 	
if [ "$1" == "ab1_clean" ]; then    ######LISTO
	ab1_clean.sh $data_path
fi

if [ "$1" == "char" ]; then   ######LISTO TODO MENOS DFAST que lo hemos ejecutado con el dfast_clean.sh includo 
	. ~soft_bio_267/initializes/init_autoflow  
	
	 rm -r $genome_analysis_path/char/dfast_0000  #$genome_analysis_path/char/pyani_0000  $genome_analysis_path/char/Sibelia_0000 $data_path/total_genomes $data_path/all_genome_list
	 mkdir -p $data_path/total_genomes 
	 mkdir -p $genome_analysis_path/char
	
	 ln -s $data_path/genomes_problem/*fasta $data_path/total_genomes
	 ln -s $data_path/genomes_down/genomes/*.fasta $data_path/total_genomes

	 ls $data_path/genomes_down/genomes > $data_path/genome_candidate_list

	 ls $data_path/genomes_problem > $data_path/all_genome_list
	 ls $data_path/genomes_down/genomes >> $data_path/all_genome_list
	
	vars=`echo "
		\\$data_path=$project_path'/data',
		\\$ab1_name=$data_path/ab1_name,
		\\$scripts_path=$script_path
		" | tr -d [:space:]`

	AutoFlow -c 1 -s -w $template_path/genome_analysis.af -V $vars -t "1-20:00:00" -o $genome_analysis_path/char
	
fi

##################################
## TP FLOW
##################################

if [ "$1" == "protein_db" ]; then  ######LISTO

	protein_db.sh  $script_path $data_path $SCRATCH/genome_analysis "Shewanella"
fi 

if [ "$1" == "tp_case" ]; then  ######LISTO

	. ~soft_bio_267/initializes/init_autoflow 
	rm -r $genome_analysis_path/transposon/executions
	output=$genome_analysis_path/transposon/executions

	while read genome 
	do 
		mkdir -p $output/$genome
		vars=`echo "
			\\$genome_name=$genome,
			\\$data_path=$project_path'/data',
			\\$genome_analysis_path=$genome_analysis_path,
			\\$genome_seq=$data_path/total_genomes/$genome,
			\\$prot_database=$data_path/tp_data
			" | tr -d [:space:]`
		AutoFlow -c 1 -s -w $template_path/tpflow -V $vars -o "$output/"$genome $2
	done < $data_path/all_genome_list
	
fi 

##################################
## TarSynFlow
##################################

if [ "$1" == "genes" ]; then

	output=$genome_analysis_path/genes_identification/Tarsynflow

	rm -r $output
	mkdir -p $output/proteome

	sbatch $sh_file_path/reduce_prot_redundancy.sh  $data_path $genome_analysis_path 

fi

if [ "$1" == "genes_comps" ]; then

	make_all_comps.sh $data_path $script_path $genome_analysis_path $template_path

fi

if [ "$1" == "genes_results_protein" ]; then

	sbatch $sh_file_path/get_all_results.sh $genome_analysis_path/genes_identification/Tarsynflow/results $data_path $genome_analysis_path/genes_identification/Tarsynflow/comps  $script_path get_protein
fi 

if [ "$1" == "genes_results_annotation" ]; then

	get_all_results.sh $genome_analysis_path/genes_identification/Tarsynflow/results $data_path $genome_analysis_path/genes_identification/Tarsynflow/comps  $script_path get_annotations
fi  

if [ "$1" == "seqs" ]; then

	extract_seqs.sh "$genome_analysis_path/genes_identification/Tarsynflow/comps/Shewanella_Pdp11.fasta/"`head -n 1 data/gen_refs`"/procompart_0000/coord_table_with_strand" $data_path/genomes_problem $genome_analysis_path/genes_identification/Tarsynflow $data_path $script_path 
	
fi 

##################################
## Genomic Islands
##################################


if [ "$1"  == "GI" ]; then 
	HTTP_API_token="b800318f-766a-e3ad-d20a-74b797db9969" 
	input=$genome_analysis_path/char/dfast_0000/genome_annotation
	output=$genome_analysis_path/genomic_island

	mkdir -p $output
	rm -r $output/Island_Viewer4_results
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
	rm -r $output
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
	grep -h -v '#' $genome_analysis_path/char/blastn_0000/blast_16S/blast_* | awk '{if (($3>98.7)&&($15=100)) {print $0}}'| tr ' ' '\t' | sed s'/1_S/\tS/g' | cut -f 1,2,3,15 | sed s'/.fasta//g'| sed s'/_1492R//g' | sed s'/_/ /g'> $results_path/blast_16

	### Complete genome comparison 
	grep -v "Q_Shewanella_Pdp11" $genome_analysis_path/char/pyani_0000/genome_pyani_anim/matrix_identity_1.tab | grep -v "R_scaffold" | cut -f 1-9 | sed s'/Shewanella / S./g' | sort >> $results_path/pyani_identity
	sed -i '1ishewanella strains \t Pdp11 \t SH12 \t SH16 \t SH4 \t SH6 \t SH9 \t SdM1 \t SdM2' $results_path/pyani_identity
	grep -v "Q_Shewanella_Pdp11" $genome_analysis_path/char/pyani_0000/genome_pyani_anim/matrix_coverage_1.tab | grep -v "R_scaffold" | cut -f 1-9 | sed s'/Shewanella / S./g' | sort >> $results_path/pyani_coverage
	sed -i '1ishewanella strains \t Pdp11 \t SH12 \t SH16 \t SH4 \t SH6 \t SH9 \t SdM1 \t SdM2' $results_path/pyani_coverage 
	grep -v "shewanella strains" $results_path/pyani_identity | sed s'/ /_/g' | awk '{if (($2>0.9)||($3>0.9)||($4>0.9)||($5>0.9)||($6>0.9)||($7>0.9)||($8>0.9)||($9>0.9)) {print $0}}' | sed s'/_/ /g' | cut -f 1 > $results_path/pyani_identity_top
	
	###### Synteny block and Sibelia 
	cp $genome_analysis_path/char/Sibelia_0000/Shewanella_Pdp11.fasta/Shewanella_baltica_128.fasta/circos/circos.png $results_path/S_baltica_128:Pdp11.png 
	cp $genome_analysis_path/char/Sibelia_0000/Shewanella_Pdp11.fasta/Shewanella_putrefaciens_4H.fasta/circos/circos.png $results_path/S_putrefaciens_4H:Pdp11.png

	### Genome annotation by DFAST
	cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/*/COG_percent_tab |  sed s'/\tS/\nS/g' | sed s'/\tun/\nun/g' | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/COG_percent
	sed -i '1ishewanella strains \t %' $results_path/COG_percent
	cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/Total_cog_table | sed s'/Shewanella_/S./g'| sed s'/_/ /g' > $results_path/COG_annotation
	cat $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/Total_cog_table_relative | sed s'/Shewanella_/S./g'| sed s'/_/ /g' > $results_path/COG_annotation_relative
	
	### Transposon 
	awk '{if ($2>0) {print $0}}' $genome_analysis_path/transposon/executions/Total_absolute | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/Total_absolute
	sed -i '1ishewanella strains \t transposable elements' $results_path/Total_absolute  
	awk '{if ($2>0) {print $0}}' $genome_analysis_path/transposon/executions/Total_relative | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/Total_relative
	sed -i '1ishewanella strains \t transposable elements' $results_path/Total_relative

	less -S $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,2,3,4 | tr ',' '\t' | cut -f 1,2,3,4 > $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/Pdp11_tp_1
	less -S $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/summary.txt | cut -f 1,5 | tr ',' '\t' | cut -f 1,2 > $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/Pdp11_tp_2
	merge_tabular.rb $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/Pdp11_tp_1 $genome_analysis_path/transposon/executions/Shewanella_Pdp11.fasta/transposons_finder.rb_0000/results/Pdp11_tp_2 | sed s'/Q_//g' | sed s'/_/ /' > $results_path/Pdp11_tp
	
	### Transposon matrix 
	cat $genome_analysis_path/transposon/executions/Tab_interrupt | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/Tp_interrupt
	cat $genome_analysis_path/transposon/executions/Tab_interrupt_number | sed s'/Shewanella_/S./g' | sed s'/_/ /g'>  $results_path/Tp_interrupt_total
	cat $genome_analysis_path/transposon/executions/Tab_transposase | sed s'/Shewanella_/S./g' | sed s'/_/ /g'>  $results_path/Tp_transposable
	cat $genome_analysis_path/transposon/executions/Tab_transposase_number  | sed s'/Shewanella_/S./g' | sed s'/_/ /g'> $results_path/Tp_transposable_total

	### Tarsynflow
	grep -v 'Entry' $genome_analysis_path/genes_identification/Tarsynflow/results/Ref_specific_filtered_annotated_prots | cut -f 1,5,7 > $results_path/specific_genes
	
	### Genomic Island Total
	cat $genome_analysis_path/genomic_island/genomic_island_results/Total_GI | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/GI_total
	sed -i '1ishewanella strains \t GIs number' $results_path/GI_total

	### Genomic Island COG categories
	grep "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_name
	grep -v "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_value
	merge_tabular.rb  $data_path/COG_categories_all $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_value | cut -f 2-139 > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI
	cat $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_name $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/enrichment_GI_cat

	### Genomic Island COG categories
	grep "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_name
	grep -v "category" $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_value
	merge_tabular.rb  $data_path/COG_categories_all $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_value | cut -f 2-139 > $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_relative
	cat $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_category_relative_name $genome_analysis_path/genomic_island/genomic_island_results/enrichment_GI_relative | sed s'/Shewanella_/S./g' | sed s'/_/ /g' > $results_path/enrichment_GI_cat_relative


	### Prophage 
	#cat $genome_analysis_path/prophage/Total_phage_top| sed s'/Shewanella/S./g' > $results_path/Total_phage
	#sed -i '1ishewanella strains \t prophages number' $results_path/Total_phage

	########################################
	############REPORT
	########################################

	#$results_path/Total_phage
	paths=`echo -e "
	$results_path/blast_16,
	$results_path/COG_percent,
	$results_path/COG_annotation,
	$results_path/COG_annotation_relative,
	$results_path/pyani_identity,
	$results_path/pyani_coverage,
	$results_path/Total_absolute,
	$results_path/Total_relative,
	$results_path/Pdp11_tp,
	$results_path/Tp_interrupt,
	$results_path/Tp_transposable,
	$results_path/specific_genes,
	$results_path/GI_total,
	$results_path/enrichment_GI_cat,
	$results_path/enrichment_GI_cat_relative
	" | tr -d [:space:]`
	report_html -t $template_path/report_template.erb -d $paths -o $results_path/project_report
	
fi
