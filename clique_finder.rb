require_relative "exon.rb"
require_relative "organism.rb"
require_relative "clique.rb"

class CliqueFinder

	attr_accessor :organisms
	attr_accessor :cliques

	def initialize(organisms, file_path = "")
		self.organisms = organisms
		self.cliques = []
	end

	def find_cliques
		clique_number = 0
		clique_data = "organism_name, exon_index,cliques\n"
		self.organisms.each do |organism|
			organism.exons.each_with_index do |exon, exon_index|
				next if !exon.cliques.empty?
				exon_uids = get_exon_uids_for_clique(exon)
				exons_in_clique = mark_exons_in_clique( exon, exon_uids, clique_number )
				cliques << Clique.new(clique_number, exons_in_clique)
				clique_number += 1
			end
		end
		self.organisms.each do |organism|
			organism.exons.each_with_index do |exon, exon_index|
				clique_data += "#{organism.name},#{exon_index+1},#{exon.cliques.join(' ')}\n"
			end
		end
		File.open("cliques.txt", "w") { |file| file.write(clique_data) }
	end

private

	def get_exon_uids_for_clique(exon)
		exon_connected_uids = get_connected_uids(exon)
		exon.connections.each do |connected_exon|
			exon_connected_uids = exon_connected_uids & get_connected_uids(connected_exon)
		end
		return exon_connected_uids
	end

	def mark_exons_in_clique(exon, uids_in_clique, clique_number)
		marked_exons = []
		exon.cliques << clique_number
		exon.connections.each do |connected_exon|
			if uids_in_clique.include?(connected_exon.uuid)
				connected_exon.cliques += [clique_number]
				marked_exons << connected_exon
			end
		end
		return marked_exons
	end

	def get_connected_uids(exon)
		return (exon.connections.map(&:uuid) + [exon.uuid])
	end



end