require "matrix"
require_relative "cat_cat_matcher.rb"
require_relative "data_processor.rb"

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

class LocalResult

  attr_accessor :start_positions
  attr_accessor :end_positions
  attr_accessor :align_1
  attr_accessor :align_2
  attr_accessor :score

  def initialize(results)
    self.start_positions = results["start_positions"]
    self.end_positions = results["end_positions"]
    self.align_1 = results["align_1"]
    self.align_2 = results["align_2"]
    self.score = results["score"]
  end

  def from_start_1?
    return start_1 == 0
  end

  def from_start_2?
    return start_2 == 0
  end

  def till_end_1?(seq)
    return end_1 == (seq.length-1)
  end

  def till_end_2?(seq)
    return end_2 == (seq.length-1)    
  end

  def start_1
    return self.start_positions[0]
  end

  def start_2
    return self.start_positions[1]
  end

  def end_1
    return self.end_positions[0]
  end

  def end_2
    return self.end_positions[1]
  end

  def aligns
    return [align_1, align_2]
  end

end

class CatCatLocalAligner

  attr_accessor :seq_1
  attr_accessor :seq_2

  attr_accessor :blosum


  def initialize(seq_1 = "", seq_2 = "", blosum = {})
    self.seq_1 = seq_1
    self.seq_2 = seq_2
    self.blosum = blosum.empty? ? DataProcessor.parse_blosum : blosum
  end

  def align
    alignement_matrix = Matrix.zero(self.seq_1.length+1,self.seq_2.length+1)
    alignement_matrix, back_ways = fill_matrix(alignement_matrix)
    return LocalResult.new(make_back_way(alignement_matrix, back_ways))
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
        left = matrix[i,j-1] + blosum["A"]["-"]
        top = matrix[i-1,j] + blosum["A"]["-"]
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
    end_positions = [max_position[0]-1, max_position[1]-1]
    while true
      break if matrix[max_position[0], max_position[1]] == 0
      current_way = ways[max_position.join('_')]
      aligned_1 += current_way[0] == 0 ? "-" : self.seq_1[max_position[0]-1]
      aligned_2 += current_way[1] == 0 ? "-" : self.seq_2[max_position[1]-1]
      max_position = [max_position,current_way].transpose.map {|x| x.reduce(:+)}
    end
    start_positions = max_position
    align_1 = aligned_1.reverse
    align_2 = aligned_2.reverse
    results = { 
      "start_positions" => start_positions,
      "end_positions" => end_positions,
      "align_1" => align_1,
      "align_2" => align_2,
      "score" => count_score(align_1, align_2)
      }
    return results
  end

  def count_score(str_1, str_2)
    result = 0
    start_gap = true
    str_1.chars.each_with_index do |char, index|
      if char == "-" || str_2[index] == "-"
        if start_gap
          blosum_value = self.blosum["A"]["-"]
          start_gap = false
        else
          blosum_value = -2
        end 
      else

        blosum_value = self.blosum[char]
        start_gap = true
      end
      if !blosum_value.nil?
        result += self.blosum[char][str_2[index]] 
      else
        result += self.blosum["A"]["-"]
      end

    end
    return result
  end

end