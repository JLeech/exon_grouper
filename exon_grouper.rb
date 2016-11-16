require "csv"
require "json"

require_relative "exon.rb"
require_relative "organism.rb"
require_relative "data_processor"
require_relative "exon_matcher.rb"
require_relative "group_saver.rb"
require_relative "clique_finder.rb"

class ExonGrouper
  
  attr_accessor :path_to_file
  attr_accessor :path_to_allignement
  attr_accessor :percent
  attr_accessor :organism_number

  attr_accessor :blossum_matrix

  attr_accessor :organisms
  attr_accessor :max_group
  attr_accessor :output_filename
  attr_accessor :exons_hash

  def initialize(options = {})

    self.path_to_file = options["file"]
    self.path_to_allignement = options["allignement"]
    self.percent = options["percent"].to_i
    self.organism_number = options["organism_number"].to_i
    self.output_filename = options["output_filename"]

    self.blossum_matrix = DataProcessor.parse_blossum
    self.organisms = []
    self.exons_hash = {}
  end

  def prepare_data
    # отбираются только первые self.organism_number организмов
    all_organisms = DataProcessor.new(self.path_to_file, self.path_to_allignement).prepare
    self.organisms = all_organisms[0..(self.organism_number-1)]
    #self.organisms = [all_organisms[0],all_organisms[2],all_organisms[3],all_organisms[4],all_organisms[8],all_organisms[9]]
    clear_output_file
  end

  def group
    time1 = Time.now
    exons = []
    # создаются связи между экзонами (какие вложены в какие)
    make_connections
    # reallocate_connections
    self.organisms.each do |organism|
      exons += organism.exons
    end
    time2 = Time.now
    puts "making connections: #{time2 - time1}"
    #нумеруются группы вложенности
    make_groups(exons)
    time3 = Time.now
    puts "making groups: #{time3 - time2}"
    
  end

  # def make_cliques
  #     clique_finder = CliqueFinder.new(self.organisms, output_filename)
  #     clique_finder.find_cliques
  #     clique_finder.print_cliques
  # end

  def make_connections
    ExonMatcher.csv_header(output_filename)
    ExonMatcher.clear_output_file(output_filename)
    pair_counter = 0
    self.organisms.each_with_index do |organism, index|
      puts "#{index} : #{organism.name}"
      break if index+1 == self.organisms.length
      self.organisms[(index+1)..-1].each do |match_organism|
        organism.exons.each_with_index do |exon, organism_exon_index|
          exon_includes = false # флаг того, вложен ли экзон в кого-то
          match_organism.exons.each_with_index do |match_exon, match_organism_index|
            # проверяем вложен ли экзон, с учётом процента совпадения из options
            if exon.include?(match_exon, self.percent)
              sequences = [exon.allignement, match_exon.allignement]
              coords = [] << exon.get_coords << match_exon.get_coords
              sequences_data = {pair_id: get_pair_id(pair_counter),
                        org_name: organism.name,
                        exon_index: exon.uuid,
                        match_org_name: match_organism.name,
                        match_exon_index: match_exon.uuid,
                        }
              exon_matcher = ExonMatcher.new(sequences, coords, sequences_data, blossum_matrix,
                               organism, match_organism, exon, match_exon)
              exon_matcher.count_everything
              borders = get_borders(exon_matcher, exon, match_exon, sequences_data)
              File.open("#{self.output_filename}_borders.csv", 'a') { |file| file.write(borders) }
              if ([exon_matcher.rloc_1, exon_matcher.rloc_2].max > 0.2)
                exon_includes = true
                exon.connections << match_exon
                match_exon.connections << exon
                set_exon_coefs(exon, match_exon, exon_matcher)
                exon.local_borders << [exon_matcher.local_data['start_position_1'], exon_matcher.local_data['end_position_1']]
                match_exon.local_borders << [exon_matcher.local_data['start_position_2'], exon_matcher.local_data['end_position_2']]
                File.open("#{self.output_filename}_graph_borders.csv", 'a') { |file| file.write(borders) }
              end
              exon_matcher.print_for_csv(output_filename)
              exon_matcher.print_statistics_for_txt(output_filename)
              pair_counter += 1
            
            elsif exon_includes
              break
            end
          end
        end
      end
    end
    puts "pairs: #{pair_counter}"
  end

  def make_groups(exons)
    exons.each do |exon|
      self.exons_hash[exon.uuid] = exon  
    end
    
    group = group_green(exons)
  end

  def group_green(exons)
    group = 0
    exons.each do |exon|
      if ((exon.group.empty?) & (exon.green?))
        cf = CliqueFinder.new(exon)
        cf.make_graph
        cliques = cf.make_cliques
        group = mark_cliques(cliques, group)
      end
    end
    return group
  end

  def mark_cliques(cliques, group)
    cliques.each do |clique|
      clique.each do |exon_in_clique_id|
        self.exons_hash[exon_in_clique_id].group << group
      end
      group += 1
    end
    return group
  end

  def set_exon_coefs(exon, match_exon, exon_matcher)
    rloc_max = [exon_matcher.rloc_1, exon_matcher.rloc_2].max
    local_length_coef_max = [ exon_matcher.local_data["local_length_1_coef"], exon_matcher.local_data["local_length_2_coef"] ].max
    min_local_length = [ exon_matcher.local_data["raw_length_1"], exon_matcher.local_data["raw_length_2"] ].min
    matching_letters_length_coef = [exon_matcher.matching_letters/exon.alignement_no_gap_length,exon_matcher.matching_letters/match_exon.alignement_no_gap_length].max
    exon.min_local_lengths[match_exon.organism_index] << min_local_length
    exon.r_maxes[match_exon.organism_index] << rloc_max
    exon.local_length_coef_maxes[match_exon.organism_index] << local_length_coef_max
    exon.matching_letters[match_exon.organism_index] << matching_letters_length_coef

    match_exon.min_local_lengths[exon.organism_index] << min_local_length
    match_exon.r_maxes[exon.organism_index] << rloc_max
    match_exon.local_length_coef_maxes[exon.organism_index] << local_length_coef_max
    match_exon.matching_letters[exon.organism_index] << matching_letters_length_coef
  end

  def draw_as_svg_rectangels(data_to_show)
    svg_width = self.organisms.first.allignement_length*2 + 200
    svg_height = self.organisms.count * 40 + 40
    output_file_name = self.output_filename + "_#{data_to_show}"
    File.open("#{output_file_name}.svg", 'w') do |file|
      file.write("<svg width=\"#{svg_width+100}\" height=\"#{svg_height}\">")
      file.write("<rect x=\"0\" y=\"0\" width=\"#{svg_width+100}\" height=\"#{svg_height}\" style=\"fill:white;\" />")
      self.organisms.each_with_index do |organism, index|
        draw_organism_line(index, svg_width, file)
        print_organism_name(organism, index, file)
        organism.exons.each do |exon|
          draw_exon_box(index, exon, file, exon.send(data_to_show))
        end
      end
      file.write("</svg>")
    end
    puts "exon number : #{organisms.map{ |org| org.exons.length }.inject(:+)}"
    `inkscape -z -e #{output_file_name}.png -w #{svg_width} -h #{svg_height} #{output_file_name}.svg`
  end

  def print_groups_to_csv
    group_saver = GroupSaver.new(self.organisms, output_filename)
    group_saver.save_to_csv
  end
  
  def get_borders(aligner, exon, match_exon, sequence_data)
    data1 = [exon.uuid, match_exon.uuid, sequence_data[:pair_id],aligner.local_data['start_position_1'],aligner.local_data['end_position_1'],aligner.local_data['start_position_2'],aligner.local_data['end_position_2']]
    data1 = data1.join(",") + "\n"
    data2 = [match_exon.uuid, exon.uuid, sequence_data[:pair_id],aligner.local_data['start_position_2'],aligner.local_data['end_position_2'],aligner.local_data['start_position_1'],aligner.local_data['end_position_1']]
    data2 = data2.join(",") + "\n"
    return data1+data2
  end

private

  def clear_output_file
    header = "ex1,ex2,pair_id,st1,end1,st2,end2\n"
    File.open("#{self.output_filename}_borders.csv", 'w') { |file| file.write(header) }
    header = "ex1,ex2,pair_id,st1,end1,st2,end2\n"
    File.open("#{self.output_filename}_graph_borders.csv", 'w') { |file| file.write(header) }
  end

  def get_pair_id(number)
    return "A"+"0"*(5-number.to_s.length) + number.to_s
  end

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

  def draw_exon_box(index, exon, file, data_to_show)
    y_coords = 40*(index+1)
    x_start_coords = 100
    width = exon.finish - exon.start
    color = exon.get_svg_color
    file.write("<rect x=\"#{(exon.start+x_start_coords)*2}\" y=\"#{y_coords-15}\" width=\"#{width*2}\" height=\"30\" style=\"fill:#{color};stroke:black;stroke-width:1\" />\n")
    #file.write("<text x=\"#{(exon.start+x_start_coords)*2+10}\" y=\"#{y_coords}\" fill=\"black\">(#{exon.start}:#{exon.finish})</text>")
    file.write("<text x=\"#{(exon.start+x_start_coords)*2}\" y=\"#{y_coords}\" fill=\"black\">#{data_to_show}</text>")
  end

end