class Exon

	attr_accessor :start
	attr_accessor :finish
	attr_accessor :allignement

	attr_accessor :connections
#	attr_accessor :real_connections
	attr_accessor :group
#	attr_accessor :cliques

	attr_accessor :r_maxes
	attr_accessor :local_length_coef_maxes
	attr_accessor :min_local_lengths

	attr_accessor :organism_index
	attr_accessor :exon_index

	attr_accessor :uuid

	attr_accessor :local_borders

	def initialize(start, finish, allignement, organism_index = 1, exon_index = 1)
		self.start = start
		self.finish = finish
		self.allignement = allignement
		self.connections = []
#		self.real_connections = []
		self.group = -1
		self.uuid = (organism_index+1)*100 + exon_index + 1
		# self.cliques = []
		self.organism_index = organism_index
		self.exon_index = exon_index
		self.local_borders = []
		self.r_maxes = []
		self.local_length_coef_maxes = []
		self.min_local_lengths = []
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

	# def collect_cliques
	# 	if !real_connections.empty?
	# 		self.cliques = real_connections.map(&:cliques).flatten.uniq
	# 	end
	# end

	def max_blossum(blossum)
		max_score = 0.0
		allignement.split("").each { |char| max_score += blossum[char][char] }
		return max_score
	end

	# def has_local_overlap?
	# 	return false if local_borders.empty?
	# 	starts = self.local_borders.map(&:first)
	# 	ends = self.local_borders.map(&:last)
	# 	if starts.sort.last > ends.sort.first
	# 		return true
	# 	else
	# 		return false
	# 	end
	# end

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
		puts allignement
		puts "-----------------"
	end

	def get_svg_color
		color = "rgb(18, 221, 232)"
		color_letter = rmax_leng_coef
        color = "rgb(115, 219, 67)" if color_letter == "g"
        color = "rgb(232, 221, 18)" if color_letter == "y"
        return color
	end

	def get_validations
		counter_min = 0.0
		counter_max = 0.0
		self.r_maxes.each_with_index do |r_max, index|
			if ((r_max > 0.3) & (self.local_length_coef_maxes[index] > 0.75) & (self.min_local_lengths[index] > 5) )
				counter_max += 1.0
			elsif ((r_max < 0.3) & (self.min_local_lengths[index] < 5) )
				counter_min += 1.0
			end
		end
		return [counter_min, counter_max]
	end

private

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

end