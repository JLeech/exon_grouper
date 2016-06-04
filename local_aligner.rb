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
		self.seq_1 = seq_1
		self.seq_2 = seq_2
		self.blosum = LocalAligner.parse_blosum if blosum.empty?
	end

	def align
		alignement_matrix = Matrix.zero(self.seq_1.length+1,self.seq_2.length+1)
		alignement_matrix, back_ways = fill_matrix(alignement_matrix)
		return make_back_way(alignement_matrix, back_ways)
	end

	def fill_matrix(matrix)
		back_ways = {}
		ways = [[-1,-1],[0,-1],[-1,0],[0,0]]
		for i in 1..(matrix.row_count-1)
			for j in 1..(matrix.column_count-1)
				diagonal = matrix[i-1,j-1]+blosum[seq_1[i-1]][seq_2[j-1]]
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

	def count_score(str_1, str_2)
		result = 0
		str_1.chars.each_with_index do |char, index|
			result += self.blosum[char][str_2[index]]
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

seq_1 = "PPDGAPQDVQLEAISSQGIKVTWK"
seq_2 = "SPDGPPQEVQLEALSSQSVKVTWK"

la = LocalAligner.new(seq_1, seq_2)
puts seq_1.length
puts la.align
#x = Matrix[[1, 10, 5], [3, 4,0]]
#puts "#{x.get_max_position}"

