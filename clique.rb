require_relative "exon.rb"

class Clique
	
	attr_accessor :clique_number
	attr_accessor :exons


	def initialize(clique_number, exons)
		self.clique_number = clique_number
		self.exons = exons
	end

	def length
		return exons.length
	end

end