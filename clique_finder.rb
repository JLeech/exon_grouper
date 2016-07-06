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
		time_start = Time.now
		puts "clique finder"
		graph = construct_graph
		graph_constructed = Time.now
		puts "graph constructed #{graph_constructed - time_start}"
		self.cliques = kerbosh([], graph.keys, [], graph)
		cliques_found = Time.now
		puts "cliques found #{cliques_found - graph_constructed}"
		mark_exons_in_clique(cliques)
		exons_marked = Time.now
		puts "exons marked #{exons_marked - cliques_found}"
		#collect_cliques
	end

	def collect_cliques
		self.organisms.each { |organism| organism.exons.map(&:collect_cliques) }
	end

	def print_cliques
		
		organism_exon_clique = ["organism,org_index, exon_index, exon_id, cliques\n"]
		self.organisms.each_with_index do |organism, org_index|
			organism.exons.each_with_index do |exon, exon_index|
				data = "#{organism.name},#{org_index+1},#{exon_index+1},#{exon.uuid},#{exon.cliques.join(',')}"
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

end