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
		# margin allignement

		current_allignement, match_allignement = margin_allignemet(current_allignement, match_allignement, [start, finish], [exon.start, exon.finish], blossum)
		# puts "xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
		puts "#{current_allignement}"
		puts "#{match_allignement}"
		# puts "xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
		start_gap_flag = true

		current_allignement.split("").each_with_index do |current_char, position|
			match_char = match_allignement[position]
			if current_char == "-" || match_char == "-"
				match_score += start_gap_flag ? -9 : -2
				start_gap_flag = false
			else
				start_gap_flag = true
				match_score += blossum[current_char][match_char].nil? ? -9 : blossum[current_char][match_char]
			end

		end
		return match_score

	end

	def print
		puts "-----------------"
		puts "#{@start} - #{@finish}"
		puts "#{@group}"
		puts allignement
		puts "-----------------"
	end

private

	def margin_allignemet(seq_1, seq_2, coords_1, coords_2, blossum)
		start_diff = coords_1[0] - coords_2[0]
		if start_diff > 0
			seq_2 = "-"*start_diff + seq_2
		else
			seq_1 = "-"*(-start_diff) + seq_1
		end

		end_diff = coords_1[1] - coords_2[1]
		if end_diff > 0
			seq_2 = seq_2 + "-"*((end_diff - start_diff).abs)
		else
			seq_1 = seq_1 + "-"*((end_diff - start_diff).abs)
		end
		return seq_1, seq_2
	end



end