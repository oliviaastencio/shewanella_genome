#!/usr/bin/env ruby
require 'optparse'
def do_makeblastdb(seqs, output, dbtype)
	cmd="makeblastdb -in - -out #{output} -max_file_sz '100000GB' -title #{File.basename(output)} -dbtype #{dbtype} -parse_seqids"
	puts cmd
	IO.popen(cmd,'w+') {|makedb|
		makedb.sync = TRUE
		makedb.write(seqs)
		makedb.close_write
		puts makedb.readlines
		makedb.close_read
	}
end

options = {}

optparse = OptionParser.new do |opts|
  options[:name] = nil
  opts.on( '-n', '--name STRING', 'Database name in case the creation of a local DB') do |dname|
	options[:name] = dname
  end

  options[:user_fasta] = nil
  opts.on( '-f', '--user_fasta FILE', 'Use a custom fasta file to build the user database') do |file|
		options[:user_fasta] = file
  end

  # Set a banner, displayed at the top of the help screen.
  opts.banner = "Usage: #{File.basename(__FILE__)} [options]  \n\n"

  # This displays the help screen
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end

end # End opts

# parse options and remove from ARGV
optparse.parse!
#PATHS
user_db_folder = File.join(Dir.pwd, options[:name])
Dir.mkdir(user_db_folder) if !File.exists?(user_db_folder)
output_file_path = File.join(user_db_folder, options[:name])
seqs = File.open(options[:user_fasta]).read
#MAKE
do_makeblastdb(seqs, output_file_path, 'prot')
