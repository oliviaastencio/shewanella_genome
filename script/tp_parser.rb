#! /usr/bin/env ruby

#####################################################
## METHODS
######################################################
input = ARGV[0]
output = ARGV[1]

## cheching if not nil 
if input.nil? || output.nil?
  puts "Nil: ruby script.rb input.tsv output.tsv"
  exit 1
end

## input files
table = File.readlines(input).map { |line| line.chomp.split("\t") }

if table.empty?
  puts "El archivo está vacío."
  exit 1
end


col_to_delete = []
(1...table[0].size).each do |col_idx|  
  all_zeros = table[1..-1].all? { |row| row[col_idx].to_f == 0 }  
  col_to_delete << col_idx if all_zeros
end

## write output
table.each { |row| col_to_delete.reverse_each { |col_idx| row.delete_at(col_idx) } }


File.open(output, "w") do |file|
  table.each { |row| file.puts row.join("\t") }
end

