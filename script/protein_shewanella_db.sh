#! /usr/bin/env bash
script_path=$1
data_path=$2
transposon_analysis_path=$3
keyword=$4
export PATH=$script_path:$PATH
. ~soft_cvi_114/initializes/init_fln

mkdir -p $data_path/tp_data  $transposon_analysis_path/transposon/executions

curl  'https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query='$keyword > $data_path/tp_data/total_prots.fasta  #'https://www.uniprot.org:443/uniprot/?query=%20taxonomy:'$taxonomy'&format=fasta'

current=`pwd`
cd $data_path/tp_data
which make_user_db.rb
make_user_db.rb -l -f total_prots.fasta -n local_database
cd $current

while read line ; do
	genome_seq=$data_path/genomes_problem/$line.fasta
	echo $genome_seq >> $transposon_analysis_path/transposon/executions/fna.list
done < $data_path/genome_name