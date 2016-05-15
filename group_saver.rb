require_relative "exon.rb"
require_relative "organism.rb"

class GroupSaver

	attr_accessor :organisms
	attr_accessor :exons_pointers
	attr_accessor :organisms_strings
	attr_accessor :output_filename


	def initialize (organisms, output_filename)
		self.organisms = organisms
		self.exons_pointers = Array.new(organisms.length) { |i| i = 0 }
		self.organisms_strings = Array.new(organisms.length) { |i| i = "" }
		self.output_filename = output_filename
	end

	def save_to_csv
		# неоптимально
		while true	
			group_number = get_next_group_number
			break if group_number < 0
			exon_per_group = Array.new(self.organisms.length) { |i| i = 0 }
			self.organisms.each_with_index do |organism, index|
				exon_per_group[index] = organism.exons.select { |exon| exon.group == group_number }.length
			end
			max_elems_in_group = exon_per_group.max
			exon_per_group.each_with_index do |exon_number, organism_index|
				self.organisms_strings[organism_index] += "#{group_number},"*exon_number + ","*(max_elems_in_group - exon_number) + ","
				self.exons_pointers[organism_index] += exon_number
			end
		end
		File.write("#{output_filename}_groups.csv", self.organisms_strings.join("\n") )
	end


private

	def get_next_group_number
		left_exons = []
		self.exons_pointers.each_with_index do |exon_position, index|
			exon = self.organisms[index].exons[exon_position]
			next if exon.nil?
			left_exons << exon
		end
		return -1 if left_exons.compact.empty?
		exon_coords = left_exons.map(&:start)
		exon_groups = left_exons.map(&:group)
		#exon_mean = left_exons.map(&:finish).inject(:+).to_f / left_exons.size
		left_exon_group = left_exons[ exon_coords.index(exon_coords.min) ].group
		
		not_left_exons_for_group = left_exons.select { |exon| exon.group != left_exon_group }
		left_exons_for_group = left_exons.select { |exon| exon.group == left_exon_group }
		exon_mean = left_exons_for_group.map(&:finish).inject(:+).to_f / left_exons_for_group.size
		out_group = left_exon_group
		not_left_exons_for_group.each do |exon|
			if exon.finish < exon_mean
				out_group = exon.group
				break
			end
		end
		return out_group
	end

end
