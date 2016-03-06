require "csv"
require_relative "exon.rb"
require_relative "gene.rb"
require_relative "connector.rb"
require_relative "statistic.rb"

class DataProcessor

	attr_accessor :raw_data
	attr_accessor :genes

	def initialize(raw_data)
		@raw_data = raw_data
		@genes = []
	end

	def prepare
		@raw_data.each_with_index do |row,index|
			@genes <<  parse_gene(row,index)
		end
		return @genes
	end

	def parse_gene(row,index)
		organism_name = row.first.split(";").first.to_s
		current_start = row.first.split(";").last.gsub!("(","").to_i
		row.delete_at(0)
		current_exons = []
		row.each_with_index do |exon_range, exon_index|
			exon_finish,exon_start = exon_range.split(";")
			current_exons.push(Exon.new(current_start,exon_finish.gsub(")","").to_i, "#{index} #{exon_index}"))
			current_start = exon_start.gsub("(","").to_i unless exon_start.nil?
		end
		gene = Gene.new(organism_name,current_exons,index)
		gene.count_distances
		return gene
	end
end

class ExonVisualizer

	attr_accessor :path_to_file
	attr_accessor :raw_data
	attr_accessor :genes

	def initialize(path_to_file)
		@path_to_file = path_to_file
		@genes = []
	end

	def read_csv
		@raw_data = CSV.read(@path_to_file)
	end

	def prepare_data
		@genes = DataProcessor.new(@raw_data).prepare
		Connector.new(@genes).make_connections(85)
	end

	def make_json
		json_main = "{ \"nodes\": [ \n"
		@genes.each do |gene|
			gene.exons.each_with_index do |exon,index|
				json_main += exon.node_as_json
				json_main += ",\n" unless index == gene.exons.length
			end
		end
		json_main += "]}"
		File.write('exons.json',json_main)
	end

	def make_svg
		puts "wanna svg"
		file = 'exons.svg'
		File.write(file," ")
		width = 2200
		height = 16500
		svg_main = "<html>\n<body>\n<svg width=\"#{width}\" height=\"#{height}\" fill=\"rgb(0,0,0)\" >\n"
		svg_main += "<rect width=\"#{width}\" height=\"#{height}\" fill=\"rgb(250,250,250)\" stroke-width=\"1\" stroke=\"rgb(0,0,0)\" /> \n"
		write_part(file,svg_main)
		svg_main = draw_connections
		write_part(file,svg_main)
		svg_main = draw_nodes
		write_part(file,svg_main)
		svg_main = draw_nodes_coords
		write_part(file,svg_main)
		svg_main = draw_gene_names
		write_part(file,svg_main)
		svg_main = "</svg >\n</body>\n</html>\n"
		write_part(file,svg_main)

	end

	def write_part(file,svg_main)
		open(file, 'a') { |f|
  		f.puts svg_main
		}
	end

	def draw_nodes
		svg_main = ""
		@genes.each do |gene|
			gene.exons.each_with_index do |exon,index|
				svg_main += exon.node_as_svg
				svg_main += ",\n" unless index == gene.exons.length
			end
		end
		return svg_main
	end

	def draw_connections
		svg_main = ""
		@genes.each do |gene|
			gene.exons.each do |exon|
				svg_main += exon.connections_as_svg
			end
			svg_main += "\n"
		end
		return svg_main
	end

	def draw_nodes_coords
		svg_main = ""
		@genes.each do |gene|
			gene.exons.each do |exon|
				svg_main += exon.coords_as_svg
			end
			svg_main += "\n"
		end
		return svg_main
	end

	def draw_gene_names
		svg_main = ""
		@genes.each do |gene|
			svg_main += gene.name_as_svg
		end
		return svg_main
	end

	def statistic
		stat = Statistic.new(@genes)
		stat.intersections
		stat.intersections_percent
	end

end

exon_viz = ExonVisualizer.new("DSCAM_coords.csv")
exon_viz.read_csv
exon_viz.prepare_data
#exon_viz.make_svg
exon_viz.statistic



