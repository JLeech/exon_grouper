class Exon

	attr_accessor :start
	attr_accessor :finish
	attr_accessor :connections
	attr_accessor :group

	def initialize(start, finish, graph_number)
		@start = start
		@finish = finish
		@connections = []
		@group = -1
	end

	def include?(exon, match_persent)
		first_range = (start..finish)
		second_range = (exon.start..exon.finish)
		match_length = ([first_range.begin,second_range.begin].max..[first_range.max,second_range.max].min).size
		if match_length >= [first_range.size, second_range.size].min * match_persent / 100.0
			return true
		else
			return false
		end
	end

	def intersection_length(exon)
		first_range = (start..finish)
		second_range = (exon.start..exon.finish)
		match_length = ([first_range.begin,second_range.begin].max..[first_range.max,second_range.max].min).size
		return match_length
	end
end