require_relative "exon.rb"

class Organism

	attr_accessor :name
	attr_accessor :exons
	attr_accessor :number
	attr_accessor :allignement_length

	def initialize(name, exons, number, allignement_length)
		@name = name
		@exons = exons
		@number = number
		@allignement_length = allignement_length
		puts name
		puts allignement_length
	end

end