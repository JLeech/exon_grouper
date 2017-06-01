class Exon

	attr_accessor :start
	attr_accessor :finish
	attr_accessor :alignment

	attr_accessor :organism_index
	attr_accessor :exon_index
	attr_accessor :color_letter
	attr_accessor :uuid


	def initialize(start, finish, alignment, organism_index = 1, exon_index = 1)
		self.start = start
		self.finish = finish
		self.alignment = alignment
		self.uuid = (organism_index+1)*100 + exon_index + 1
		self.organism_index = organism_index
		self.exon_index = exon_index
	end


	def get_exons_matching_coords(exon)
		first_range = (start..finish)
		second_range = (exon.start..exon.finish)
		match_length = ([first_range.begin, second_range.begin].max..[first_range.max, second_range.max].min).size
		return match_length
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

	def get_coords_for_table
		return "(#{get_coords.join(';')})"
	end

	def as_range
		return (start..finish)
	end

	def get_index_in_org
		return self.uuid%100 - 1
	end

	def print
		puts "-----------------"
		puts "#{self.start} - #{self.finish}"
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

private



end