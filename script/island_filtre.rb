#! /usr/bin/env ruby

######METHODS#################

def clean_island(input)
	island_length = []
	File.open(input).each do |line|
		line = line.split(",")
		if line[2].to_i > 10000 
			if !line[10].include?('hypothetical protein')
				island_length << line.join("\t")
			end
		end 
	end 
	return island_length
end 

def write_island_results(output, island_length)
	File.open(output, 'w') do |f1|
		f1.puts island_length
	end 
end 

##############OPTPARSE############
input, output = ARGV   
#####MAIN#####################
island_length = clean_island(input)
write_island_results(output, island_length)