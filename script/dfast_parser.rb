#! /usr/bin/env ruby

#####################################################
## METHODS
######################################################
def save_record(genes, cog_text, tigr_text, gene_name)    
    cog_id = nil 
    cog_category = nil
    tigr_id = nil
    if cog_text != ''
        cog_id = cog_text[/COG\d+/]
        cog_category_selection = cog_text[/Category:\w+/] 
	cog_category_filtre = cog_category_selection.split(':')   
        cog_category = cog_category_filtre[1].split('')     
	end
    if tigr_text != '' 
            tigr_id = tigr_text[/TIGR\d+/] 
    end
    genes << [gene_name, cog_id, cog_category, tigr_id] if !cog_id.nil? || !tigr_id.nil? 
end

def load_dfast_file(input)
    genes = []
    gene_name = nil
    cog_text = ''
    tigr_text = ''
    processing_tigr = false
    processing_cog = false

    #UPLOAD FILES
    File.open(input).each do |line| 
        line.chomp! 
           if line.include?('/locus_tag=')
                if !gene_name.nil?
                    save_record(genes, cog_text, tigr_text, gene_name)
                    tigr_text = ''
                    cog_text = ''     
                end
                fields = line.split('=')
                gene_name = fields[1].gsub('"', '')
            elsif line.include?('/note=')
                    processing_tigr = false 
                    processing_cog = false 
                    note = line.split('=') 
                    if note[1].include?('TIGR:')
                            tigr_text << note[1]
                            processing_tigr = true 
                elsif note[1].include?('COG:')
                cog_text << note[1]
                processing_cog = true
                    end 
            elsif processing_tigr == true
                    tigr_text << line.strip
        elsif processing_cog == true
                cog_text << line.strip 
        end
    end
    save_record(genes, cog_text, tigr_text, gene_name) 
    return genes 
end

def compute_cog_statistics(genes)
    cog_table = [] 
    total_category = Hash.new(0) 
    genes.each do |gene|
        cog_categories = gene[2]
        if !cog_categories.nil? 
            cog_categories.each do |cog_category|
                total_category[cog_category] += 1
            end
        end
    end 
    category_list = total_category.to_a.sort{|a1, a2| a1[0] <=> a2[0]}
    category_list.each do |category, count|
       cog_table << [category,count]
    end 
    return cog_table
end

def write_gene_table(output, genes)
    File.open(output, 'w') do |f2|
        genes.each do |gene|
            gene[2] = gene[2].join(',') if !gene[2].nil?
            f2.puts gene.join("\t")
        end 
    end
end

def write_cog_table(cog_table_path, cog_table)
    File.open(cog_table_path, 'w') do |f1|
        cog_table.each do |cog|
            f1.puts cog.join("\t")
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
genes = load_dfast_file(input)
cog_table = compute_cog_statistics(genes)
cog_table_path = File.join(File.dirname(output), 'cog_table.txt')
write_cog_table(cog_table_path, cog_table)
write_gene_table(output, genes)