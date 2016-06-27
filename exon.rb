require "securerandom"
class Exon

	attr_accessor :start
	attr_accessor :finish
	attr_accessor :allignement

	attr_accessor :connections
	attr_accessor :real_connections
	attr_accessor :group
	attr_accessor :cliques

	attr_accessor :organism_index
	attr_accessor :exon_index

	attr_accessor :uuid

	def initialize(start, finish, allignement, organism_index = 1, exon_index = 1)
		self.start = start
		self.finish = finish
		self.allignement = allignement
		self.connections = []
		self.real_connections = []
		self.group = -1
		#self.uuid = SecureRandom.hex(10)
		self.uuid = (organism_index+1)*100 + exon_index + 1
		self.cliques = []
		self.organism_index = organism_index
		self.exon_index = exon_index
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

	def click_include(exon)
		return true if start == exon.start && finish == exon.finish
		return false
	end

	def collect_cliques
		if connections.empty? && !real_connections.empty?
			self.cliques = real_connections.map(&:cliques).flatten.uniq
		end
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

end