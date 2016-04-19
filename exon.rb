class Exon

	attr_accessor :start
	attr_accessor :finish
	attr_accessor :connections
	attr_accessor :group
	attr_accessor :allignement

	def initialize(start, finish, allignement)
		@start = start
		@finish = finish
		@allignement = allignement
		@connections = []
		@group = -1
		@click = -1
	end

	def include?(exon, match_persent, blossum)
		# процент вложенности считается для наименьшего экзона
		if match_persent == 100
			#count_with_blossum(exon, blossum)/max_blossum(blossum)
			return click_include(exon) 
		end
		first_range = (start..finish)
		second_range = (exon.start..exon.finish)
		match_length = ([first_range.begin, second_range.begin].max..[first_range.max, second_range.max].min).size
		if match_length >= [first_range.size, second_range.size].min * match_persent / 100.0
			#puts count_with_blossum(exon, blossum)/max_blossum(blossum)
			return true
		else
			return false
		end
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

	def count_with_blossum(exon, blossum)
		match_score = 0.0
		current_allignement = allignement
		match_allignement = exon.allignement
		if start < exon.start
			match_allignement = "$"*(exon.start - start) + match_allignement
		elsif start > exon.start
			current_allignement = "$"*(start - exon.start) + current_allignement
		end
		current_allignement.split("").each_with_index do |current_char, index|
			break if match_allignement[index].nil? || current_char.nil?
			break if match_allignement[index].strip.empty? || current_char.strip.empty?
			if blossum[current_char][match_allignement[index]].nil?
				match_score += 0
			else
				match_score += blossum[current_char][match_allignement[index]]
			end
		end
		return match_score

	end 

end