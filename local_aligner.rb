require "matrix"

class Matrix
	
	def []=(i,j,k)
		@rows[i][j]=k
	end # добавление метода записи в матрицы, а то они по дефолту immutable
	
	def get_max_position
		max = -9999
		position = []
		@rows.each_with_index do |row, row_index|
			row.each_with_index do |val, column_index|
				if val > max
					position = [row_index, column_index]
					max = val
				end
			end
		end
		
		return position
	end
	
end

class LocalAligner


	attr_accessor :seq_1
	attr_accessor :seq_2
	attr_accessor :blosum

	BLOSUM_PATH = "./data/blosum_penalty_matrix"

	def initialize(seq_1 = "", seq_2 = "", blosum = {})
		self.seq_1 = seq_1.gsub("-","")
		self.seq_2 = seq_2.gsub("-","")
		self.blosum = blosum.empty? ? LocalAligner.parse_blosum : blosum

	end

	def align
		alignement_matrix = Matrix.zero(self.seq_1.length+1,self.seq_2.length+1)
		alignement_matrix, back_ways = fill_matrix(alignement_matrix)
		return format_results( make_back_way(alignement_matrix, back_ways))
	end

	def fill_matrix(matrix)
		back_ways = {}
		ways = [[-1,-1],[0,-1],[-1,0],[0,0]]
		for i in 1..(matrix.row_count-1)
			for j in 1..(matrix.column_count-1)
				blosum_value = blosum[seq_1[i-1]][seq_2[j-1]]
				if blosum_value.nil?
					diagonal = matrix[i-1,j-1]	
				else
					diagonal = matrix[i-1,j-1]+blosum_value
				end
				left = matrix[i,j-1]+blosum["A"]["-"]
				top = matrix[i-1,j]+blosum["A"]["-"]
				max = [diagonal, left, top, 0].max
				back_ways["#{i}_#{j}"] = ways[[diagonal, left, top, 0].index(max)]
				matrix[i,j] = max
			end
		end
		return [matrix, back_ways]
	end

	def make_back_way(matrix, ways)
		max_position = matrix.get_max_position
		aligned_1 = ""
		aligned_2 = ""
		start_positions = [max_position,[-1,-1]].transpose.map {|x| x.reduce(:+)}
		while true
			break if matrix[max_position[0],max_position[1]] == 0
			current_way = ways[max_position.join('_')]
			aligned_1 += self.seq_1[max_position[0]-1] == 0 ? "-" : self.seq_1[max_position[0]-1]
			aligned_2 += self.seq_2[max_position[1]-1] == 0 ? "-" : self.seq_2[max_position[1]-1]
			max_position = [max_position,current_way].transpose.map {|x| x.reduce(:+)}
		end
		end_positions = max_position
		results = { 
			"start_positions" => start_positions,
			"end_positions" => end_positions,
			"align_1" => aligned_1.reverse,
			"align_2" => aligned_2.reverse,
			"score" => count_score(aligned_1, aligned_2)
		  }
		return results
	end

	def format_results(results)
		start_position_1 = results["start_positions"][0]
		start_position_2 = results["start_positions"][1]
		end_position_1 = results["end_positions"][0]
		end_position_2 = results["end_positions"][1]
		length_1 = end_position_1 - start_position_1
		length_2 = end_position_2 - start_position_2
		score_1 = count_score(results["align_1"],results["align_1"])
		score_2 = count_score(results["align_2"],results["align_2"])
		local_length_coef = ([length_1,length_2].min.to_f)/([length_1,length_2].max.to_f)
		local_score_1_coef = score_1.to_f/results["score"]
		local_score_2_coef = score_2.to_f/results["score"]
		formatted_result = {
			"start_position_1" => start_position_1,
			"start_position_2" => start_position_2,
			"end_position_1" => end_position_1,
			"end_position_2" => end_position_2,
			"score_1" => score_1,
			"score_2" => score_2,
			"local_length_coef" => local_length_coef,
			"local_score_1_coef" => local_score_1_coef,
			"local_score_2_coef" => local_score_2_coef,
			"align_1" => results["align_1"],
			"align_2" => results["align_2"]
		}
		return formatted_result
	end

	def count_score(str_1, str_2)
		result = 0
		str_1.chars.each_with_index do |char, index|
			blosum_value = self.blosum[char][str_2[index]]
			unless blosum_value.nil?
				result += self.blosum[char][str_2[index]]	
			end
		end
		return result
	end

	def print_matrix(mat)
		for i in 0..(mat.row_count-1)
			for j in 0..(mat.column_count-1)
				printf " #{mat[i,j]}"
			end
			puts
		end
	end 

	def self.parse_blosum
		alphabet = []
		tmp_array = []
		blosum_matrix = Hash.new { |hash, key| hash[key] = Hash.new { |hash, key| hash[key] = 0 } }
		File.readlines(BLOSUM_PATH).each_with_index do |line, index|
			next if [0,2].include?(index)
			if index == 1
				alphabet = line.strip.split(" ")
				next
			end

			values_array = line.strip.split(" ")
			alphabet.length.times { |pos| tmp_array += [alphabet[pos], values_array[pos].to_i] }
			blosum_matrix[alphabet[index - 3]] = Hash[*tmp_array]
			tmp_array = []

		end
		return blosum_matrix
	end

end

# seq_1 = "PPDGAPQDVQLEAISSQGIKVTWK"
# seq_2 = "SPDGPPQEVQLEALSSQSVKVTWK"

# la = LocalAligner.new(seq_1, seq_2)
# puts seq_1.length
# puts la.align
# #x = Matrix[[1, 10, 5], [3, 4,0]]
# #puts "#{x.get_max_position}"

