require "csv"
require "json"

require_relative "exon.rb"
require_relative "gene.rb"
require_relative "data_processor"

class ExonGrouper
	
	attr_accessor :path_to_file
	attr_accessor :percent
	attr_accessor :organism_number

	attr_accessor :raw_data
	attr_accessor :genes
	attr_accessor :max_group

	def initialize(options)

		@path_to_file = options["file"]
		@percent = options["percent"].to_i
		@organism_number = options["organism_number"].to_i

		@genes = []
	end

	def read_csv
		@raw_data = CSV.read(@path_to_file)
	end

	def prepare_data
		# отбираются только первые @organism_number организмов
		@genes = DataProcessor.new(@raw_data).prepare[0..(@organism_number)]
	end

	def group
		exons = []
		# создаются связи между экзонами (какие вложены в какие)
		make_connections
	  @genes.each do |gene|	
	  	exons += gene.exons
	  end
	  # нумеруются группы вложенности
	  make_groups(exons)
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

	def make_connections
		@genes.each_with_index do |gene, index|
			break if index+1 == @genes.length-1
			@genes[(index+1)..-1].each do |match_gene|
				gene.exons.each do |exon|
					exon_includes = false # флаг того, вложен ли экзон в кого-то
					match_gene.exons.each_with_index do |match_exon,current_position|
						# проверяем вложен ли экзон, с учётои процента совпадения из options
						if exon.include?(match_exon, @percent) 
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

	def print_groups_coords
		max_exon_count = @genes.map { |gene| gene.exons.length }.max
		puts max_exon_count
		@genes.each do |gene|
			exon_groups = " "
			gene.exons.each do |exon|
				exon_groups += " \"#{exon.start}:#{exon.finish}\", "
			end
			(max_exon_count - gene.exons.count).times do 
				exon_groups += " \"_____\", "
			end
			puts "[ #{exon_groups}],"
		end
	end

	def print_groups
		max_exon_count = @genes.map { |gene| gene.exons.length }.max
		puts max_exon_count
		@genes.each do |gene|
			exon_groups = ""
			gene.exons.each do |exon|
				exon_groups += "#{exon.group}, "
			end
			(max_exon_count - gene.exons.count).times do 
				exon_groups += "-1, "
			end
			puts "[#{exon_groups}],"
		end
	end

	def print_group_count
		@genes.each do |gene|
			groups = Array.new(@max_group, 0)
			gene.exons.each do |exon|
				groups[exon.group] += 1
			end
			groups_string = "[ "
			groups.each do |val|
				groups_string += "#{val}, "
			end
			groups_string += "], "
			puts groups_string
		end
	end

	def print_groups_as_csv
		CSV.open("groups.csv", "wb") do |csv|
			out_data = ["organism name \\ exons"]
			csv << out_data
		  @genes.each do |gene|
				out_data = []
				out_data += [gene.name]
				out_data += gene.exons.map(&:group)		  	
				csv << out_data	
		  end
		end
	end

	def print_group_count_as_csv
		CSV.open("groups_count.csv", "wb") do |csv|
			out_data = ["organism name"]
			@max_group.times do |iter|
				out_data += ["group #{iter}"]
			end
			csv << out_data
			@genes.each do |gene|
				groups = Array.new(@max_group, 0)
				gene.exons.each do |exon|
					groups[exon.group] += 1
				end	
				out_data = [gene.name]
				out_data += groups
				csv << out_data
			end
		end
	end

end