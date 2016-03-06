require_relative "exon.rb"

class Gene

	attr_accessor :name
	attr_accessor :exons
	attr_accessor :number

	def initialize(name, exons, number)
		@name = name
		@exons = exons
		@number = number
	end

end