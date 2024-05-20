#! /usr/bin/env ruby
def load_dfast_file(input)
    categories = nil
    category_table = [] 
    total_category = Hash.new(0) 
    File.open(input).each do |line|
    line.chomp!
    categories = line.split('')
        categories.each do |category|
                total_category[category] += 1
            end
        end 
        category_list = total_category.to_a.sort{|a1, a2| a1[0] <=> a2[0]} 
        category_list.each do |category, count|
        category_table << [category,count]
        end  
    return category_table
end

def write_category_table(cat_table_path, category_table)
    File.open(cat_table_path, 'w') do |f1|
        category_table.each do |cat|
            f1.puts cat.join("\t")
        end 
    end 
end 


#####################################################
## OPTPARSE
######################################################
input, output = ARGV                               
#####################################################
## MAIN
######################################################

category_table = load_dfast_file(input)
cat_table_path = File.join(File.dirname(output), 'category_table.txt')
write_category_table(cat_table_path, category_table)
 