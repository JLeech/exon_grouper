require_relative "exon.rb"
require_relative "exon_grouper.rb"

class ExonMatcher

	attr_accessor :blosum
	attr_accessor :seq_1
	attr_accessor :seq_2
	attr_accessor :coords_1
	attr_accessor :coords_2


	def initialize
		self.blosum = ExonGrouper.new.blossum_matrix
	end

	def parse_seq
		File.readlines("./data/for_exon_matcher").each_with_index do |line, index|
			self.seq_1 = line.strip if index == 0
			self.seq_2 = line.strip if index == 1
			self.coords_1 = line.strip.split(":").map(&:to_i) if index == 2
			self.coords_2 = line.strip.split(":").map(&:to_i) if index == 3
		end
	end

	def margin_allignements
		start_diff = coords_1[0] - coords_2[0]
		if start_diff > 0
			self.seq_2 = "-"*start_diff + self.seq_2
		else
			self.seq_1 = "-"*(-start_diff) + self.seq_1
		end

		end_diff = coords_1[1] - coords_2[1]
		if end_diff > 0
			self.seq_2 = self.seq_2 + "-"*((end_diff - start_diff).abs)
		else
			self.seq_1 = self.seq_1 + "-"*((end_diff - start_diff).abs)
		end
	end

	def count_blosum_and_print
		score_line = ""
		seq_1_line = ""
		seq_2_line = ""
		score = 0
		start_gap_flag = true
		self.seq_1.split("").each_with_index do |seq_1_char, position|
			seq_2_char = seq_2[position]
			if seq_1_char == "-" || seq_2_char == "-"
				score += start_gap_flag ? -9 : -2
				start_gap_flag = false
			else
				start_gap_flag = true
				score += blosum[seq_1_char][seq_2[position]]
			end

			score_line += "#{score.to_s}|"
			seq_1_line += " "*score.to_s.length << seq_1_char
			seq_2_line += " "*score.to_s.length << seq_2_char

		end
		puts "score: #{score}"
		puts "_______________________"
		puts score_line
		puts seq_1_line
		puts seq_2_line
	end

	def print_allignements
		puts "#{seq_1}"
		puts "#{seq_2}"
		puts "_______________________"
	end

	def print_coords
		puts "seq_1 coords: #{coords_1}"
		puts "seq_2 coords: #{coords_2}"
		puts "_______________________"
	end

end


exon_matcher = ExonMatcher.new
exon_matcher.parse_seq

exon_matcher.print_coords
exon_matcher.print_allignements
exon_matcher.margin_allignements
exon_matcher.print_allignements

exon_matcher.count_blosum_and_print


