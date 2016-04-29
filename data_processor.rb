require_relative "exon.rb"
require_relative "organism.rb"

class DataProcessor

	attr_accessor :raw_data
	attr_accessor :organisms

	def initialize(path_to_data_file, path_to_allignement_file)
		@raw_data = CSV.read(path_to_data_file)
		@path_to_allignement_file = path_to_allignement_file
		@organisms = []	
	end

	def prepare
		@raw_data.each_with_index do |row, index|
			@organisms <<  parse_organism(row, index)
		end
		return @organisms
	end

private

	def parse_organism(row, index)
		current_exons = []
		parsed_coords_and_name = row.join(",").split(";")
		coordinates = parsed_coords_and_name[1..(-1)]
		organism_name = parsed_coords_and_name.first.strip
		organism_allignement = get_allignement_for_organism(organism_name)
		# парсит csv и сохраняет экзоны в организмы
		coordinates.each do |exon_coordinates|
			exon_start, exon_finish = get_coords(exon_coordinates)
			current_exons.push( Exon.new(exon_start, exon_finish, organism_allignement[exon_start..exon_finish]) )
		end
		organism = Organism.new(organism_name, current_exons, index, organism_allignement.length)
		return organism
	end

	def get_coords(coords_block)
		return coords_block.gsub("(","").gsub(")","").gsub(" ","").split(",").map(&:to_i)
	end

	def get_allignement_for_organism(organism_name)
		found_flag = false
		allignement_string = ""
		File.readlines(@path_to_allignement_file).each do |line|
			if found_flag
				break if line.start_with?(">")
				allignement_string += line.strip
			else
				found_flag = true if line.start_with?(">#{organism_name}")
			end
		end
		return allignement_string
	end
end