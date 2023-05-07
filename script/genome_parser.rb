#! /usr/bin/env ruby

##########################################################
##METHODS##
#######################################################
def load_file_genome(input)
	genomes = []
	processing_genome = false  
	File.open(input).each do |line|
	line.chomp!
		if line.include?('complete')
			processing_genome = false
			if !line.include?('plasmid')
				genomes << line
				processing_genome = true
			end
		elsif processing_genome == true
				genomes << line.strip
		end
	end
	return genomes
end

def load_file_plasmid(input)
	plasmids = []
	processing_plasmid = false
	File.open(input).each do |line|
	line.chomp!
		if line.include?('complete')
			processing_plasmid = false
			if line.include?('plasmid')
				plasmids << line 
				processing_plasmid = true
			end 
		elsif processing_plasmid == true
			plasmids << line.strip  
		end
	end
	return plasmids
end

def write_genome_file(output, genomes)
	File.open(output, 'w') do |f1|
		genomes.each do |genome|
			f1.puts genome
		end
	end 
end 

def write_plasmid_file(output_1, plasmids)
 	File.open(output_1, 'w') do |f1|
 		plasmids.each do |plasmid|
 			f1.puts plasmid
 		end
	end 
end 

###################################################################
###OPTPARSE
###################################################################
input, output, output_1 = ARGV 

################################################################
###MAIN
###############################################################
genomes = load_file_genome(input)
plasmids = load_file_plasmid(input)
write_genome_file(output, genomes)
write_plasmid_file(output_1, plasmids)
