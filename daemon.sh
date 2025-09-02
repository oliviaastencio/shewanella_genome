#! /usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out


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

if [ "$1" == "char" ]; then   ######LISTO pendiente pyani
	. ~soft_bio_267/initializes/init_autoflow  

	 rm -r  $genome_analysis_path/char/pyani_0000 $data_path/total_genomes $data_path/all_genome_list # $genome_analysis_path/char/dfast_0000 $genome_analysis_path/char/Sibelia_0000 
	 mkdir -p $data_path/total_genomes 
	 mkdir -p $genome_analysis_path/char
	
	 ln -s $data_path/genomes_problem/*fasta $data_path/total_genomes
	 ln -s $data_path/genomes_down/genomes/*.fasta $data_path/total_genomes

	 ls $data_path/genomes_down/genomes > $data_path/genome_candidate_list

	 ls $data_path/genomes_problem > $data_path/all_genome_list
	 ls $data_path/genomes_down/genomes >> $data_path/all_genome_list

	current=`pwd`
	cd $data_path/16S_NCBI/
		module load blast_plus/2.15.0+
		makeblastdb -in 16S_ribosomal_RNA.fasta -parse_seqids -dbtype nucl -out 16S_db
	cd $current

	vars=`echo "
		\\$data_path=$project_path'/data',
		\\$ab1_name=$data_path/ab1_name,
		\\$scripts_path=$script_path
		" | tr -d [:space:]`

	AutoFlow -c 1 -s -w $template_path/genome_analysis.af -V $vars -t "2-00:00:00" -o $genome_analysis_path/char
	
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
	done < $data_path/all_genome_list #all_genome_list

fi 

if [ "$1" == "tp_matrix" ]; then

	sbatch $sh_file_path/Tp_matrix.sh $data_path $genome_analysis_path $genome_analysis_path/transposon/executions $genome_analysis_path/transposon/results

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


if [ "$1"  == "GI_up" ]; then 

	GI.sh UP $data_path  $genome_analysis_path
fi 

if [ "$1" == "GI_result" ]; then ###PENDIENTE

	GI.sh result $data_path  $genome_analysis_path
fi 

if [ "$1" == "GI_clean" ]; then ###PENDIENTE

	sbatch $sh_file_path/GI.sh clean $data_path  $genome_analysis_path  
	
fi 

##################################
## Prophage by PHASTEST
##################################
if [ "$1" == "phage_visorter" ]; then 

	sbatch $sh_file_path/phage.sh $data_path $genome_analysis_path
fi 

if [ "$1" == "results" ]; then  
	
 	results.sh $data_path $results_path $genome_analysis_path Shewanella Shewanella_putrefaciens_ATCC_8071 Shewanella_baltica_OS678_OS678
fi 
## REPORTING
####################################
# /mnt/home/users/pab_001_uma/pedro/html_reporting/template
if [ "$1" == "report" ]; then 
	. ~soft_bio_267/initializes/init_ruby
	
	########################################
	############REPORT
	########################################
	
	paths=`echo -e "
		$results_path/report_img/blast_16,
		$results_path/report_img/report_Total_cog_table, 
		$results_path/report_img/report_Total_cog_table_relative,
		$results_path/report_img/pyani_matrix_identity,
		$results_path/report_img/pyani_matrix_coverage,
		$results_path/report_img/report_Tp_absolute,
		$results_path/report_img/report_Tp_relative,
		$results_path/report_img/Tp_Pdp11,
		$results_path/report_img/report_Tp_interrupt,
		$results_path/report_img/report_Tp_transposase,
		$results_path/report_img/specific_genes,
		$results_path/report_img/report_Total_GI,
		$results_path/report_img/report_Total_GI_relative,
		$results_path/report_img/report_enrichment_GI_category,
		$results_path/report_img/report_enrichment_GI_category_relative,
		$results_path/report_img/report_Total_phage,
		$results_path/report_img/report_Total_phage_relative
		" | tr -d [:space:]`
	report_html -t $template_path/report_template.erb -d $paths -o $results_path/project_report
	
fi
