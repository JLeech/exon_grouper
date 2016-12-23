# require_relative "exon.rb"
# require_relative "organism.rb"

# class GroupSaver

# 	attr_accessor :organisms
# 	attr_accessor :output_filename


# 	def initialize (organisms, output_filename)
# 		self.organisms = organisms
# 		self.output_filename = output_filename
# 	end

# 	def save_to_csv
# 		save_org_exon
# 		save_org_clique
# 	end

# 	def save_org_exon
# 		max_exons_length = organisms.map(&:exons).map(&:length).max-1
# 		lines =  ["org name\\exon index;" + (0..max_exons_length).to_a.join(";")]
# 		organisms.each do |org|
# 			lines << org.name + ";" + org.exons.map(&:group).map{ |ex| ex.join(",") }.join(";")
# 		end
# 		File.write("#{output_filename}_org_exon_groups.csv", lines.join("\n") )
# 	end

# 	def save_org_clique
# 		max_clique_num = organisms.map(&:exons).flatten.map(&:group).flatten.max
# 		lines = []
# 		organisms.each do |org|
# 			lines << org.name + ";"
# 		end
# 		org_clique_arr = []
# 		organisms.each_with_index do |org, index|
# 			cliques = Hash.new(0)
# 			org.exons.map(&:group).flatten.each do |group|
# 				cliques[group] += 1
# 			end
# 			org_clique_arr[index] = cliques
# 		end

# 		(0..max_clique_num).to_a.each do |clique_number|
# 			org_clique_arr.each_with_index do |org, index|
# 				exons_in_clique = org[clique_number]
# 				lines[index] += "#{exons_in_clique};"
# 			end
# 		end
# 		lines = ["org name\\clique index;" + (0..max_clique_num).to_a.join(";")] + lines
# 		File.write("#{output_filename}_exons_in_cliques_count.csv", lines.join("\n") )

# 	end

# private

# end
