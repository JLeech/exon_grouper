require_relative "exon.rb"
require_relative "gene.rb"

class DataProcessor

	attr_accessor :raw_data
	attr_accessor :genes

	def initialize(raw_data)
		@raw_data = raw_data
		@genes = []
	end

	def prepare
		@raw_data.each_with_index do |row,index|
			@genes <<  parse_gene(row,index)
		end
		return @genes
	end

	def parse_gene(row,index)
		organism_name = row.first.split(";").first.to_s
		current_start = row.first.split(";").last.gsub!("(","").to_i
		row.delete_at(0)
		current_exons = []
		row.each_with_index do |exon_range, exon_index|
			exon_finish,exon_start = exon_range.split(";")
			current_exons.push(Exon.new(current_start,exon_finish.gsub(")","").to_i, "#{index} #{exon_index}"))
			current_start = exon_start.gsub("(","").to_i unless exon_start.nil?
		end
		gene = Gene.new(organism_name,current_exons,index)
		return gene
	end
end