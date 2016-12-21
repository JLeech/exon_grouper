require_relative "exon.rb"

class Organism

	attr_accessor :name
	attr_accessor :code_name
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
		self.code_name = name.split("_").map { |val| val[0..2] }.join("_")
		puts "#{name} : #{code_name}"
	end

end





