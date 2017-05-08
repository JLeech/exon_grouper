class Exon

	attr_accessor :start
	attr_accessor :finish
	attr_accessor :alignment

	attr_accessor :connections
#	attr_accessor :real_connections
	attr_accessor :group
#	attr_accessor :cliques

	attr_accessor :r_maxes
	attr_accessor :local_length_coef_maxes
	attr_accessor :min_local_lengths
	attr_accessor :matching_letters

	attr_accessor :organism_index
	attr_accessor :exon_index

	attr_accessor :uuid

	attr_accessor :local_borders

	def initialize(start, finish, alignment, organism_index = 1, exon_index = 1)
		self.start = start
		self.finish = finish
		self.alignment = alignment
		self.connections = []
#		self.real_connections = []
		self.group = []
		self.uuid = (organism_index+1)*100 + exon_index + 1
		# self.cliques = []
		self.organism_index = organism_index
		self.exon_index = exon_index
		self.local_borders = []
		self.r_maxes = Hash.new {|hsh, key| hsh[key] = [] }
		self.local_length_coef_maxes = Hash.new {|hsh, key| hsh[key] = [] }
		self.min_local_lengths = Hash.new {|hsh, key| hsh[key] = [] }
		self.matching_letters = Hash.new {|hsh, key| hsh[key] = [] }
	end

	def include?(exon, match_persent)
		# процент вложенности считается для наименьшего экзона
		if match_persent == 100
			return click_include?(exon) 
		end
		first_range = (start..finish)
		second_range = (exon.start..exon.finish)
		match_length = ([first_range.begin, second_range.begin].max..[first_range.max, second_range.max].min).size
		if match_length >= [first_range.size, second_range.size].min * match_persent / 100.0
			#return false if [first_range.size, second_range.size].max > 3*[first_range.size, second_range.size].min			
			return true
		else
			return false
		end
	end

	def connected_organisms
		self.connections.map(&:organism_index)
	end

	def get_exons_matching_coords(exon)
		first_range = (start..finish)
		second_range = (exon.start..exon.finish)
		match_length = ([first_range.begin, second_range.begin].max..[first_range.max, second_range.max].min).size
		return match_length
	end

	def click_include?(exon)
		return true if start == exon.start && finish == exon.finish
		return false
	end

	def get_connected_hash
		connects = Hash.new {|hsh, key| hsh[key] = [] }
		connections.each do |connected|
			connects[connected.organism_index] << connected
		end
		return connects
	end

	def max_blossum(blossum)
		max_score = 0.0
		alignment.split("").each { |char| max_score += blossum[char][char] }
		return max_score
	end

	def alignement_no_gap_length
		counter = 0.0
		self.alignment.each_char do |char|
			counter += 1 if char != "-"
		end
		return counter
	end

	def get_coords
		return [start, finish]
	end

	def get_index_in_org
		return self.uuid%100 - 1
	end

	def print
		puts "-----------------"
		puts "#{self.start} - #{self.finish}"
		puts "#{self.group}"
		puts alignment
		puts "-----------------"
	end

	def get_svg_color
		color = "rgb(18, 221, 232)"
        color = "rgb(115, 219, 67)" if green?
        color = "rgb(232, 221, 18)" if yellow?
        return color
	end

	def green?
		return color_letter == "g"
	end

	def blue?
		return color_letter == "b"
	end

	def yellow?
		return color_letter == "y"
	end

	def get_validations
		counter_min = 0.0
		counter_max = 0.0
		# con_org_id connected_organism_id
		self.r_maxes.keys.each do |con_org_id|
			if (self.r_maxes[con_org_id].max > 0.3) & (self.local_length_coef_maxes[con_org_id].max > 0.75) & (self.min_local_lengths[con_org_id].min > 5)
				counter_max += 1.0
			elsif ((self.r_maxes[con_org_id].max < 0.3) & (self.min_local_lengths[con_org_id].min < 5) )
				counter_min += 1.0
			end
		end 

		return [counter_max, counter_min]
	end

	def rmax_leng_coef
 		counter_max, counter_min = get_validations
		if counter_max/self.r_maxes.length > 0.75
			return "g"
		end
		if counter_min/self.r_maxes.length > 0.75
			return "y"
		end
		return "b"
	end

	def color_letter
		return @color_letter unless @color_letter.nil?
		@color_letter = rmax_leng_coef
		return @color_letter
	end

	def get_unique_groups
		return self.group.uniq
	end

	def not_grouped?
		return group.empty?
	end

private



end