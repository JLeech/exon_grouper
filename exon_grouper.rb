require "csv"
require "json"

require_relative "exon.rb"
require_relative "organism.rb"
require_relative "data_processor"
require_relative "exon_matcher.rb"

class ExonGrouper
    
    attr_accessor :path_to_file
    attr_accessor :path_to_allignement
    attr_accessor :percent
    attr_accessor :organism_number

    attr_accessor :blossum_matrix

    attr_accessor :organisms
    attr_accessor :max_group
    attr_accessor :output_filename

    def initialize(options = {})

        @path_to_file = options["file"]
        @path_to_allignement = options["allignement"]
        @percent = options["percent"].to_i
        @organism_number = options["organism_number"].to_i
        @output_filename = options["output_filename"]

        @blossum_matrix = DataProcessor.parse_blossum
        @organisms = []
    end

    def prepare_data
        # отбираются только первые @organism_number организмов
        @organisms = DataProcessor.new(@path_to_file, @path_to_allignement).prepare[0..(@organism_number)]
    end

    def group
        exons = []
        # создаются связи между экзонами (какие вложены в какие)
        make_connections
        @organisms.each do |organism|
            exons += organism.exons
        end
        #нумеруются группы вложенности
        make_groups(exons)
    end

    def make_connections

    	ExonMatcher.csv_header(output_filename)
    	ExonMatcher.clear_output_file(output_filename)
 
      @organisms.each_with_index do |organism, index|
        break if index+1 == @organisms.length-1
        @organisms[(index+1)..-1].each do |match_organism|
          organism.exons.each_with_index do |exon, organism_index|
            exon_includes = false # флаг того, вложен ли экзон в кого-то
            match_organism.exons.each_with_index do |match_exon, match_organism_index|
              # проверяем вложен ли экзон, с учётом процента совпадения из options
              if exon.include?(match_exon, @percent, blossum_matrix)
                pair_id = (organism.name + organism_index.to_s + match_organism.name + match_organism_index.to_s).hash
                
                sequences = [exon.allignement, match_exon.allignement]
                coords = [] << exon.get_coords << match_exon.get_coords
                sequences_data = {pair_id: pair_id,
                									org_name: organism.name,
                									exon_index: organism_index,
                									match_org_name: match_organism.name,
                									match_exon_index: match_organism_index,
                									}
                exon_matcher = ExonMatcher.new(sequences, coords, sequences_data, blossum_matrix)
                exon_matcher.count_everything
                additional_data = {}
                exon_matcher.print_for_csv(output_filename)
                exon_matcher.print_statistics_for_txt(output_filename)

                exon_includes = true
                exon.connections << match_exon
                match_exon.connections << exon
              elsif exon_includes
                break
              end
            end
          end
        end
      end
    end

    def make_groups(exons)
        group = 0
        exons.each do |exon|
            if exon.group == -1
                exon.group = group
                group += 1
                # проходим по всем вложенным, и проставляем им группу текущего
                set_groups_for_connected(exon)
            end
        end
        @max_group = group  
    end

    def set_groups_for_connected(exon)
        exon.connections.each do |connected_exon|
            if connected_exon.group == -1
                connected_exon.group = exon.group
                set_groups_for_connected(connected_exon)
            end
        end
    end

    def draw_as_svg_rectangels
        svg_width = @organisms.first.allignement_length*2 + 200
        svg_height = @organisms.count * 40 + 40
        output_file_name = path_to_allignement.split('/').last.split("_").first
        File.open("#{output_file_name}.svg", 'w') do |file|
            file.write("<svg width=\"#{svg_width+100}\" height=\"#{svg_height}\">")
            file.write("<rect x=\"0\" y=\"0\" width=\"#{svg_width+100}\" height=\"#{svg_height}\" style=\"fill:white;\" />")
            @organisms.each_with_index do |organism, index|
                draw_organism_line(index, svg_width, file)
                print_organism_name(organism, index, file)
                organism.exons.each do |exon|
                    draw_exon_box(index, exon, file)
                end
            end
            file.write("</svg>")
        end
        #`inkscape -z -e #{output_file_name}.png -w #{svg_width} -h #{svg_height} #{output_file_name}.svg`
    end


private
    
    def print_organism_name(organism, index, file)
        y_coords = 40*(index+1) - 10
        file.write("<text x=\"10\" y=\"#{y_coords}\" fill=\"black\" font-size=\"20px\">#{organism.name}</text>")
    end

    def draw_organism_line(index, svg_width, file)
        y_coords = 40*(index+1)
        x_start_coords = 100
        x_end_coords = svg_width - 100
        file.write("<line x1=\"#{x_start_coords}\" y1=\"#{y_coords}\" x2=\"#{x_end_coords}\" y2=\"#{y_coords}\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />")
    end

    def draw_exon_box(index, exon, file)
        y_coords = 40*(index+1)
        x_start_coords = 100
        width = exon.finish - exon.start
        file.write("<rect x=\"#{(exon.start+x_start_coords)*2}\" y=\"#{y_coords-15}\" width=\"#{width*2}\" height=\"30\" style=\"fill:rgb(204,255,51);stroke:black;stroke-width:1\" />\n")
        #file.write("<text x=\"#{(exon.start+x_start_coords)*2+10}\" y=\"#{y_coords}\" fill=\"black\">(#{exon.start}:#{exon.finish})</text>")
        file.write("<text x=\"#{(exon.start+x_start_coords)*2}\" y=\"#{y_coords}\" fill=\"black\">(#{exon.group})</text>")
    end

end