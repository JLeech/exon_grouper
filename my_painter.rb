require_relative "exon_grouper"

def read_options
	options = {}
	options_lines = File.read('options')
	options_lines.split("\n").each do |option_line|
		next if option_line.strip.start_with?("#")
		line_data = option_line.strip.split("=").map(&:strip)
		next if line_data.first.nil?
		options[line_data.first] = line_data.last
	end
	return options
end

def check_options(options)
	
	if options["file"].nil?
		puts "no file selected"
		return false
	end

	if options["allignement"].nil?
		puts "no allignement selected"
		return false
	end
	
	if options["percent"].nil?
		puts "no percent setted"
		return false
	else
		if options["percent"].to_i > 100 || options["percent"].to_i < 0
			puts "percent should be in range [0,100]"
			return false
		end
	end

	if options["output"].nil? || options["output"] != "csv"
		puts "format not setted or wrong. possible formats are: csv"
		return false
	end

	if options["organism_number"].nil?
		puts "more them zero organisms should be selected"
		return false
	end

	return true

end

options = read_options

exit unless check_options(options)

exon_grouper = ExonGrouper.new(options)
# создание генов и заполнение их экзонами
exon_grouper.prepare_data
# групировка на основе вложенности
exon_grouper.group
#exon_grouper.print_groups_coords
#exon_grouper.print_groups
#exon_grouper.print_group_count
if options["output"] == "csv"	
	#exon_grouper.print_group_count_as_csv
end
exon_grouper.draw_as_svg_rectangels
#exon_grouper.print_groups_coords