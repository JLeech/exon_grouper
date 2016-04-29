require "csv"
require "json"

require_relative "exon.rb"
require_relative "organism.rb"
require_relative "data_processor"

class ExonGrouper

	BLOSSUM_PATH = "./data/blosum_penalty_matrix"
	
	attr_accessor :path_to_file
	attr_accessor :path_to_allignement
	attr_accessor :percent
	attr_accessor :organism_number

	attr_accessor :blossum_matrix

	attr_accessor :organisms
	attr_accessor :max_group

	def initialize(options = {})

		@path_to_file = options["file"]
		@path_to_allignement = options["allignement"]
		@percent = options["percent"].to_i
		@organism_number = options["organism_number"].to_i

		@blossum_matrix = parse_blossum
		@organisms = []
	end

	def prepare_data
		# отбираются только первые @organism_number организмов
		@organisms = DataProcessor.new(@path_to_file, @path_to_allignement).prepare[0..(@organism_number)]
	end

	def group
		exons = []
		# создаются связи между экзонами (какие вложены в какие)
		make_connections
	  @organisms.each do |organism|	
	  	exons += organism.exons
	  end
	  #нумеруются группы вложенности
	  make_groups(exons)
	end

	def make_connections
		@organisms.each_with_index do |organism, index|
			break if index+1 == @organisms.length-1
			@organisms[(index+1)..-1].each do |match_organism|
				organism.exons.each_with_index do |exon, organism_index|
					exon_includes = false # флаг того, вложен ли экзон в кого-то
					match_organism.exons.each_with_index do |match_exon, match_organism_index|
						# проверяем вложен ли экзон, с учётом процента совпадения из options
						if exon.include?(match_exon, @percent, blossum_matrix)
							puts "_________________________"
							blos = exon.count_with_blossum(match_exon, blossum_matrix)
							max_blos = exon.max_blossum(blossum_matrix)
							puts "\n#{organism.name} : exon_number:[#{organism_index + 1}]"
							puts "#{match_organism.name} : exon_number:[#{match_organism_index + 1}]"
							puts "#{exon.start} : #{exon.finish}"
							puts "#{match_exon.start} : #{match_exon.finish}"
							puts "one to one : #{blos}"
							puts "max for one: #{max_blos}"
							puts "percents   : #{(blos/max_blos*100).round}%"
							puts "_________________________"
							exon_includes = true
							exon.connections << match_exon
							match_exon.connections << exon
						elsif exon_includes
							break
						end
					end
				end
			end
		end
	end

	def make_groups(exons)
		group = 0
		exons.each do |exon|
			if exon.group == -1
				exon.group = group
				group += 1
				# проходим по всем вложенным, и проставляем им группу текущего
				set_groups_for_connected(exon)
			end
		end
		@max_group = group  
	end

	def set_groups_for_connected(exon)
		exon.connections.each do |connected_exon|
			if connected_exon.group == -1
				connected_exon.group = exon.group
				set_groups_for_connected(connected_exon)
			end
		end
	end

	def parse_blossum
		alphabet = []
		tmp_array = []
		blossum_matrix = Hash.new { |hash, key| hash[key] = Hash.new { |hash, key| hash[key] = 0 } }
		File.readlines(BLOSSUM_PATH).each_with_index do |line, index|
			next if [0,2].include?(index)
			if index == 1
				alphabet = line.strip.split(" ")
				next
			end

			values_array = line.strip.split(" ")
			alphabet.length.times { |pos| tmp_array += [alphabet[pos], values_array[pos].to_i] }
			blossum_matrix[alphabet[index - 3]] = Hash[*tmp_array]
			tmp_array = []

		end
		return blossum_matrix
	end

	def print_groups_coords
		max_exon_count = @organisms.map { |organism| organism.exons.length }.max
		puts max_exon_count
		@organisms.each do |organism|
			exon_groups = " "
			organism.exons.each do |exon|
				exon_groups += " \"#{exon.start}:#{exon.finish}\", "
			end
			(max_exon_count - organism.exons.count).times do 
				exon_groups += " \"_____\", "
			end
			puts "[ #{exon_groups}],"
		end
	end

	def print_groups
		max_exon_count = @organisms.map { |organism| organism.exons.length }.max
		puts max_exon_count
		@organisms.each do |organism|
			exon_groups = ""
			organism.exons.each do |exon|
				exon_groups += "#{exon.group}, "
			end
			(max_exon_count - organism.exons.count).times do 
				exon_groups += "-1, "
			end
			puts "[#{exon_groups}],"
		end
	end

	def draw_as_svg_rectangels
		svg_width = @organisms.first.allignement_length*2 + 200
		svg_height = @organisms.count * 40 + 40
		output_file_name = path_to_allignement.split('/').last.split("_").first
		File.open("#{output_file_name}.svg", 'w') do |file|
			file.write("<svg width=\"#{svg_width+100}\" height=\"#{svg_height}\">")
			file.write("<rect x=\"0\" y=\"0\" width=\"#{svg_width+100}\" height=\"#{svg_height}\" style=\"fill:white;\" />")
			@organisms.each_with_index do |organism, index|
				draw_organism_line(index, svg_width, file)
				print_organism_name(organism, index, file)
				organism.exons.each do |exon|
					draw_exon_box(index, exon, file)
				end
			end
			file.write("</svg>")
		end
		#`inkscape -z -e #{output_file_name}.png -w #{svg_width} -h #{svg_height} #{output_file_name}.svg`
	end


private
	
	def print_organism_name(organism, index, file)
		y_coords = 40*(index+1) - 10
		file.write("<text x=\"10\" y=\"#{y_coords}\" fill=\"black\" font-size=\"20px\">#{organism.name}</text>")
	end

	def draw_organism_line(index, svg_width, file)
		y_coords = 40*(index+1)
		x_start_coords = 100
		x_end_coords = svg_width - 100
		file.write("<line x1=\"#{x_start_coords}\" y1=\"#{y_coords}\" x2=\"#{x_end_coords}\" y2=\"#{y_coords}\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />")
	end

	def draw_exon_box(index, exon, file)
		y_coords = 40*(index+1)
		x_start_coords = 100
		width = exon.finish - exon.start
		file.write("<rect x=\"#{(exon.start+x_start_coords)*2}\" y=\"#{y_coords-15}\" width=\"#{width*2}\" height=\"30\" style=\"fill:rgb(204,255,51);stroke:black;stroke-width:1\" />\n")
		#file.write("<text x=\"#{(exon.start+x_start_coords)*2+10}\" y=\"#{y_coords}\" fill=\"black\">(#{exon.start}:#{exon.finish})</text>")
		file.write("<text x=\"#{(exon.start+x_start_coords)*2}\" y=\"#{y_coords}\" fill=\"black\">(#{exon.finish - exon.start})</text>")
	end

	# def print_group_count
	# 	@organisms.each do |organism|
	# 		groups = Array.new(@max_group, 0)
	# 		organism.exons.each do |exon|
	# 			groups[exon.group] += 1
	# 		end
	# 		groups_string = "[ "
	# 		groups.each do |val|
	# 			groups_string += "#{val}, "
	# 		end
	# 		groups_string += "], "
	# 		puts groups_string
	# 	end
	# end

	# def print_groups_as_csv
	# 	CSV.open("groups.csv", "wb") do |csv|
	# 		out_data = ["organism name \\ exons"]
	# 		csv << out_data
	# 	  @organisms.each do |organism|
	# 			out_data = []
	# 			out_data += [organism.name]
	# 			out_data += organism.exons.map(&:group)		  	
	# 			csv << out_data	
	# 	  end
	# 	end
	# end

	# def print_group_count_as_csv
	# 	CSV.open("groups_count.csv", "wb") do |csv|
	# 		out_data = ["organism name"]
	# 		@max_group.times do |iter|
	# 			out_data += ["group #{iter}"]
	# 		end
	# 		csv << out_data
	# 		@organisms.each do |organism|
	# 			groups = Array.new(@max_group, 0)
	# 			organism.exons.each do |exon|
	# 				groups[exon.group] += 1
	# 			end	
	# 			out_data = [organism.name]
	# 			out_data += groups
	# 			csv << out_data
	# 		end
	# 	end
	# end

end