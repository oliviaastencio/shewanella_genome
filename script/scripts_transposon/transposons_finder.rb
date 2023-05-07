#!/usr/bin/env ruby


require 'optparse'
require 'report_html'

###########################################################################################################################################
## METHODS
############################################################################################################################################

def show_queries(queries)
	queries.each do |query_id, subjects|
		puts query_id
		subjects.each do |subject_id, matches|
			puts "\t#{subject_id}"
			matches.each do |match|
				puts "\t\t#{match.inspect}"
			end		
		end
	end
end

#Load blast file, select identity, query start, query stop, protein start and protein stop. We consider possible order inversion of this parameter.
def load_blast(blast_file, min_ident)
	queries = {}
	subject_lengths = {}
	File.open(blast_file).each do |line|
	
		fields = line.chomp.split("\t")
		slen = fields[13].to_f
		seq_ident = fields[2].to_f
		#seq_cover = fields[3].to_f / slen 

		next if seq_ident < min_ident #|| seq_cover < min_cover 		#Blast filtered

		query_name = fields.shift
		subject_name = fields.shift
		subject_lengths[subject_name] = slen

		q_start, q_stop, s_start, s_stop = fields[4..7].map!{|e| e.to_i} #pick hsp coordinates

		if q_start < q_stop #Direct hsp
			fields = [fields[0].to_f, q_start, q_stop, s_start, s_stop]
		elsif q_start > q_stop #reverse hsp
			fields = [fields[0].to_f, q_stop, q_start, s_start, s_stop]
		end

		check_query = queries[query_name]
		if check_query.nil?
			queries[query_name] = {subject_name => [fields] }
		else
			check_subject = check_query[subject_name]
			if check_subject.nil?
				check_query[subject_name] = [fields]
			else
				check_subject << fields
			end
		end
	end	
	return queries, subject_lengths
end

def clean_hsps(queries, subject_lengths, min_slen_coverage)
	queries2remove = []
	queries.each do |query_name, subjects|
		subjects2remove = []
		subjects.each do |subject_name, hsps|
			slen_coverage = get_hit_slen(hsps)/subject_lengths[subject_name]
			subjects2remove << subject_name if slen_coverage < min_slen_coverage || slen_coverage > 1.10 #This means that the sequences has a 10% of aa above the real length and this is an artifact
		end
		subjects2remove.each do |s_name|
			subjects.delete(s_name)
		end
		queries2remove << query_name if subjects.length == 0
	end
	queries2remove.each do |q_name|
		queries.delete(q_name)
	end
end

def get_hit_slen(hsps)
	return hsps.map{|h| h[4] - h[3]}.inject(0){|sum, i| sum + i }
end


#Divide all data: interrupted protein (next def) and putative transposon (only one match). 
def split_subjects(queries, transposon_min_length_nt)
	interrupted_proteins = {}
	putative_transposons = {}
	queries.each do |query_id, subjects|
		filtered_subject = get_interrupted_proteins(subjects, transposon_min_length_nt)
		if !filtered_subject.empty? 
			interrupted_proteins[query_id] = filtered_subject
			filtered_subject = subjects.select{|sub_id, mt| mt.length == 1}
			filtered_subject.select!{|sub_id, mt|
				put_tp = mt.first
				(put_tp[1] - put_tp[2]).abs >= transposon_min_length_nt
			}
			if !filtered_subject.empty?
				putative_transposons[query_id] = filtered_subject
			end
		end		
	end
	return interrupted_proteins, putative_transposons
end

#Process interrupted proteins: more than one match, enough gap in protein coordinates for a transposon, 
#the matches for the same protein must not overlap, the matches can have a maximum gap of 10 aa.
def get_interrupted_proteins(subjects, transposon_min_length_nt)
	filtered_subject = {}
	subjects.each do |sub_id, matches|
		if matches.length > 1									# Select matches with more than one hsp => interrupted proteins
			mt_nt_sort = matches.sort{|m1, m2| m1[1] <=> m2[1]}	# Sort by genomic coordinates
			dist_nt_between_mt = mt_nt_sort[1][1] - mt_nt_sort[0][2]						
			if dist_nt_between_mt >= transposon_min_length_nt	# Check that gap in genomic coordinates has enough room for a transposon
				mt_aa_sort = matches.sort{|m1, m2| m1[3] <=> m2[3]}	# Sort by protein coordinates
				dist_aa_between_mt = mt_aa_sort[1][3] - mt_aa_sort[0][4] # Check aa matches correlation coordinates
				if dist_aa_between_mt < 0 						# There is overlap in the protein matches
					protein_length = mt_aa_sort.map{|m| [m[3], m[4]]}.flatten.max
					if -dist_aa_between_mt.fdiv(protein_length) <= 0.1 # Matches must not overlap in protein coordinates
						filtered_subject[sub_id] = matches
					end
				elsif dist_aa_between_mt <= 10 # When matches contain a gap, this gap must have a maximum of 10 aa
					filtered_subject[sub_id] = matches
				end
			end
		end
	end
	return filtered_subject
end

#Build a new hash whose key will be the query_id and contain an array of arrays with the coordinates of the gaps of the interrupted proteins.
def detected_proteins_gaps(interrupted_proteins)
	gaps_coordinates = {}

	interrupted_proteins.each do |query_id, subjects|
		coordinates = []
		subjects.each do |subject_id, matches|
			sorted_coordinates = matches.map{|m| [m[1], m[2]]}.flatten.sort
			coordinates << [sorted_coordinates[1], sorted_coordinates[2], [subject_id]]
		end
		gaps_coordinates[query_id] = coordinates
	end		
	return gaps_coordinates
end

#Define the gaps of interrupted proteins, considering that if: 
	#Most than one proteins hace the same gap, we count this gap once.
	#Different gaps overlap, with the region overlap function, we select the biggest gap.
def define_tranposon_windows(gaps_coordinates)
	transposon_windows = {}
	gaps_coordinates.each do |query_id, gaps|
		gaps.sort!{|g1, g2| g1[1] <=> g2[1]}
		gaps = gaps.uniq
		if gaps.length > 1
			new_windows = []
			while !gaps.empty?
				ref_gap = gaps.shift
				r_start, r_stop, subject_ids = ref_gap
				gaps.each_with_index do |gap, i|
					g_start, g_stop, g_subject_ids = gap
					if region_overlap?(r_start, r_stop, g_start, g_stop)
						ref_gap = [[r_start, g_start].min, [r_stop, g_stop].max, subject_ids | g_subject_ids]
						gaps[i] = nil
					end
				end
				gaps.compact!
				new_windows << ref_gap
			end
		else
			new_windows = gaps
		end
		transposon_windows[query_id] = new_windows
	end
	return transposon_windows
end

#The putative transposons selected must be inside the gap of the interrupted proteins.
def get_transposons(putative_transposons, gaps_coordinates)
	definitive_transposons = []
	putative_transposons.each do |query_id_trans, pt_transposons|
		gaps_coordinates.each do |query_id_gaps, gaps|
			next if query_id_trans != query_id_gaps
			pt_transposons.each do |subject_id, matches|
				transposon_start = matches[0][1]
				transposon_stop = matches[0][2]
				gaps.each do |gap_start, gap_stop, interrupted_proteins|
					if region_overlap?(gap_start, gap_stop, transposon_start, transposon_stop)
						definitive_transposons << [subject_id, query_id_gaps, transposon_start, transposon_stop, [gap_start, gap_stop], interrupted_proteins]  	
					end
				end	
			end
		end
	end
	return definitive_transposons
end

#This function allow:
 #Select the biggest gap of interrupted protein, in case of the gaps overlap.
 #Select the putative transposons which are inside in the gaps of the interrupted proteins.
def region_overlap?(ref_start, ref_stop, start, stop)
	overlap = FALSE
	if (ref_start < start && ref_stop > stop) ||
		(ref_start < start && ref_stop > start) ||
		(ref_start < stop && ref_stop > stop) ||
		(ref_start >= start && ref_stop <= stop) 
		overlap = TRUE  	
	end
	return overlap
end

#Generate the reports .html.
def generate_html_report(output_paths, definitive_transposons, interrupted_proteins)
	data = []
	last_query = nil
	definitive_transposons.each do |subject_id, query_id_gaps, transposon_start, transposon_stop|
		if !last_query.nil? && query_id_gaps != last_query
			build_report(output_paths, data, interrupted_proteins, last_query)
			data = []
		else
			data << [subject_id, transposon_start, transposon_stop - transposon_start ]
		end
		last_query = query_id_gaps
	end
	build_report(output_paths, data, interrupted_proteins, last_query) if !data.empty?
end

def generate_summary_file(summary_file, definitive_transposons)
	last_chromosome = nil
	temp = []
	chromosome_matches = []
	definitive_transposons.each do |transposon_id, chromosome, transposon_start, transposon_stop, gap_coordinates, interrupted_proteins|
		if !last_chromosome.nil? && last_chromosome != chromosome
			chromosome_matches << temp
			temp = []
		end
		gap_start, gap_stop = gap_coordinates
		temp << [chromosome, gap_start, gap_stop, transposon_id, interrupted_proteins]
		last_chromosome = chromosome
	end
	chromosome_matches << temp

	File.open(summary_file, 'w') do |f|
		chromosome_matches.each do |chromosome_cluster|
			while !chromosome_cluster.empty?
				ref_trans = chromosome_cluster.first
				common_trans = chromosome_cluster.select{|trn| trn[1] == ref_trans[1] && trn[2] == ref_trans[2]}
				if !common_trans.empty?
					f.puts [ref_trans[0], ref_trans[1], ref_trans[2], common_trans.map{|t| t[3]}.join(','), 
															common_trans.map{|t| t[4]}.flatten.uniq.join(',')].join("\t")
				end
				chromosome_cluster.select!{|trn| trn[1] != ref_trans[1] || trn[2] != ref_trans[2]}
			end
		end
	end
end

def build_report(output_paths, data, interrupted_proteins, query_id)
	template = "	
		<%=
		stacked(id: :data, header: false, height: '1000px', config: {
		                'smpLabelScaleFontFactor' => 0.25
		        }
		)
		%>
	"
	subjects = interrupted_proteins[query_id]
	subjects.each do |subject_id, matches|
		matches.each do |match|
			start = match[1]
			stop = match[2]
			data << [subject_id, start, stop - start]
		end		
	end
	data.sort!{|m1, m2| m1[1] <=> m2[1]}
	container = { data: data}
	report = Report_html.new(container)
	report.build(template)
	file_path = File.join(output_paths[:html_folder], query_id + '.html')
	file_path.gsub!(/[^a-zA-Z0-9\._\/]/, '_')
	report.write(file_path)
end

###########################################################################################################################################
## OPTPARSE		
###########################################################################################################################################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_blast] = nil
  opts.on("-b", "--input_blast PATH", "Path to blast file to search tranposons") do |data| 
    options[:input_blast] = data
  end

  options[:min_ident] = 0
  opts.on("-i", "--min_ident FLOAT", "Minimum blast identity to take into account an alignment") do |data| 
    options[:min_ident] = data.to_f
  end

  options[:min_cover] = 0
  opts.on("-c", "--min_cover FLOAT", "Minimum protein coverage to take into account an alignment") do |data| 
    options[:min_cover] = data.to_f/100
  end

  options[:output_folder] = 'results'
  opts.on("-o", "--output_folder PATH", "Define path in which the program will be create a folder with the execution results") do |data| 
    options[:output_folder] = data
  end

  options[:transposon_min_length_nt] = 100
  opts.on("-t", "--transposon_min_length_nt INTEGER", "Minimum length in nucleotides of a putative tranposon to be taken into account") do |data| 
    options[:transposon_min_length_nt] = data.to_i
  end

  options[:verbose] = false
  opts.on("-v", "--verbose", "Run verbosely") do 
    options[:verbose] = true
  end

end.parse!

###########################################################################################################################################
## MAIN		
###########################################################################################################################################
output_paths = {}
output_paths[:main_folder] = options[:output_folder]
output_paths[:html_folder] = File.join(output_paths[:main_folder], 'html_reports')
output_paths[:summary_file] = File.join(output_paths[:main_folder], 'summary.txt')
output_paths[:interrupted_genes_file] = File.join(output_paths[:main_folder], 'interrupted_genes.txt')
output_paths[:definitive_transposons] = File.join(output_paths[:main_folder], 'detected_transposons.txt')

Dir.mkdir(output_paths[:main_folder]) if !File.exists?(output_paths[:main_folder])
Dir.mkdir(output_paths[:html_folder]) if !File.exists?(output_paths[:html_folder])

queries, subject_lengths = load_blast(options[:input_blast], options[:min_ident] )
clean_hsps(queries, subject_lengths, options[:min_cover])

if options[:verbose]
	puts 'ALL DATA'
	show_queries(queries)
end

interrupted_proteins, putative_transposons = split_subjects(queries, options[:transposon_min_length_nt])
if options[:verbose]
	puts 'INTERRUPTED PROTEINS'
	show_queries(interrupted_proteins)
	puts 'PUTATIVE TRANSPOSONS'
	show_queries(putative_transposons)
end
gaps_coordinates = detected_proteins_gaps(interrupted_proteins)
transposon_windows = define_tranposon_windows(gaps_coordinates)
if options[:verbose]
	puts 'GAPS COORDINATES'
	puts gaps_coordinates.inspect
end

definitive_transposons = get_transposons(putative_transposons, transposon_windows)
File.open(output_paths[:definitive_transposons], 'w') do |f|
	definitive_transposons.each do |tr_name, chr, start, stop, interrupted_prots|
		f.puts [tr_name, chr, start, stop, interrupted_prots.join(',')].join("\t")
	end
end

File.open(output_paths[:interrupted_genes_file], 'w') do |f|
	interrupted_proteins.each do |query_id, subjects|
		subjects.each do |subject_id, matches|
			transposons = definitive_transposons.select{|def_tr| def_tr.last.include?(subject_id)}.map{|m| m.first}
			next if transposons.empty?
			matches = matches.sort{|m1, m2| m1[1] <=> m2[1]}
			matches.each do |match|
				f.puts [query_id, subject_id, match[1..4], transposons.join(',')].join("\t")
			end		
		end
	end
end
generate_summary_file(output_paths[:summary_file], definitive_transposons)
generate_html_report(output_paths, definitive_transposons, interrupted_proteins)
