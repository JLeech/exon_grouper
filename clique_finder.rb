require_relative "exon.rb"
require_relative "organism.rb"
require_relative "clique.rb"

class CliqueFinder

	attr_accessor :organisms
	attr_accessor :file_path
	attr_accessor :cliques
	
	def initialize(organisms, file_path = "")
		self.organisms = organisms
		self.file_path = file_path
		self.cliques = []
	end

	def find_cliques

		graph = construct_graph
		self.cliques = kerbosh([], graph.keys, [], graph)
		mark_exons_in_clique(cliques)
		

		# clique_number = 0
		# clique_data = "organism_name, exon_index,cliques\n"
		# self.organisms.each do |organism|
		# 	organism.exons.each_with_index do |exon, exon_index|
		# 		next if !exon.cliques.empty?
		# 		exon_uids = get_exon_uids_for_clique(exon)
		# 		exons_in_clique = mark_exons_in_clique( exon, exon_uids, clique_number )
		# 		cliques << Clique.new(clique_number, exons_in_clique)
		# 		clique_number += 1
		# 	end
		# end
		# self.organisms.each do |organism|
		# 	organism.exons.each_with_index do |exon, exon_index|
		# 		clique_data += "#{organism.name},#{exon_index+1},#{exon.cliques.join(' ')}\n"
		# 	end
		# end
		# File.open("cliques.txt", "w") { |file| file.write(clique_data) }
	end

	def print_cliques
		organism_exon_clique = ["organism,org_index, exon_index, exon_id, cliques\n"]
		self.organisms.each_with_index do |organism, org_index|
			organism.exons.each_with_index do |exon, exon_index|
				data = "#{organism.name},#{org_index+1},#{exon_index+1},#{(org_index+1)*100+exon_index},#{exon.cliques.join(',')}"
				organism_exon_clique << data
			end
		end
		File.open("#{self.file_path}_org_cliques.csv", "w") { |file| file.write(organism_exon_clique.join("\n")) }
		
		clique_data = ["clique_id,exons"]
		cliques.each_with_index do |clique, clique_index|
			clique_data << "С#{clique_index}, #{clique.join(',')}"
		end
		File.open("#{self.file_path}_cliques.csv", "w") { |file| file.write(clique_data.join("\n")) }


	end

private
	
	def construct_graph
		graph = {}
		self.organisms.each do |organism|
			organism.exons.each { |exon| graph[exon.uuid] = exon.connections.map(&:uuid) }
		end
		return graph
	end

	def kerbosh(current, candidates, excluded, graph, result = [])
		result << current if candidates.empty? && excluded.empty?
		#sorted_candidates = candidates.sort { |a,b| graph[a].length <=> graph[b].length }
		current_candidates = candidates
		candidates.each do |vertex|
			result += kerbosh(current + [vertex], current_candidates & graph[vertex], excluded & graph[vertex], graph)
			current_candidates = current_candidates[1..-1]
			excluded << vertex
		end
		return result
	end

	def mark_exons_in_clique(cliques)
		exons_hash = get_exons_hash
		cliques.each_with_index do |clique, index|
			clique.each do |exon_uuid|
				exons_hash[exon_uuid].cliques += [index]
			end
		end
	end

	def get_exons_hash
		graph = {}
		self.organisms.each do |organism|
			organism.exons.each { |exon| graph[exon.uuid] = exon }
		end
		return graph
	end

	# def get_exon_uids_for_clique(exon)
	# 	exon_connected_uids = get_connected_uids(exon)
	# 	exon.connections.each do |connected_exon|
	# 		exon_connected_uids = exon_connected_uids & get_connected_uids(connected_exon)
	# 	end
	# 	return exon_connected_uids
	# end

	# def mark_exons_in_clique(exon, uids_in_clique, clique_number)
	# 	marked_exons = []
	# 	exon.cliques << clique_number
	# 	exon.connections.each do |connected_exon|
	# 		if uids_in_clique.include?(connected_exon.uuid)
	# 			connected_exon.cliques += [clique_number]
	# 			marked_exons << connected_exon
	# 		end
	# 	end
	# 	return marked_exons
	# end

	# def get_connected_uids(exon)
	# 	return (exon.connections.map(&:uuid) + [exon.uuid])
	# end

end