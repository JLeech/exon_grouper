require_relative "exon.rb"

class Organism

	attr_accessor :name
	attr_accessor :exons
	attr_accessor :number
	attr_accessor :allignement_length
	attr_accessor :allignement

	def initialize(name, exons, number, allignement)
		self.name = name
		self.exons = exons
		self.number = number
		self.allignement_length = allignement.length
		self.allignement = allignement
		puts name
	end

end