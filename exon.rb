class Exon

	attr_accessor :start
	attr_accessor :finish
	attr_accessor :connections
	attr_accessor :group
	attr_accessor :allignement
	attr_accessor :click

	def initialize(start, finish, allignement)
		self.start = start
		self.finish = finish
		self.allignement = allignement
		self.connections = []
		self.group = -1
		self.click = -1
	end

	def include?(exon, match_persent, blossum)
		# процент вложенности считается для наименьшего экзона
		if match_persent == 100
			return click_include(exon) 
		end
		first_range = (start..finish)
		second_range = (exon.start..exon.finish)
		match_length = ([first_range.begin, second_range.begin].max..[first_range.max, second_range.max].min).size
		if match_length >= [first_range.size, second_range.size].min * match_persent / 100.0
			return false if [first_range.size, second_range.size].max > 3*[first_range.size, second_range.size].min
			return true
		else
			return false
		end
	end

	def get_exons_matching_coords(exon)
		first_range = (start..finish)
		second_range = (exon.start..exon.finish)
		match_length = ([first_range.begin, second_range.begin].max..[first_range.max, second_range.max].min).size
		return match_length
	end

	def click_include(exon)
		return true if start == exon.start && finish == exon.finish
		return false
	end

	def max_blossum(blossum)
		max_score = 0.0
		allignement.split("").each { |char| max_score += blossum[char][char] }
		return max_score
	end

	def get_coords
		return [start, finish]
	end

	def print
		puts "-----------------"
		puts "#{self.start} - #{self.finish}"
		puts "#{self.group}"
		puts allignement
		puts "-----------------"
	end

private

end