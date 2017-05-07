require "csv"
require "json"
require "set"

require_relative "exon.rb"
require_relative "organism.rb"
require_relative "data_processor"
require_relative "cat_cat_matcher.rb"
require_relative "con_cat_matcher.rb"
require_relative "common.rb"

class ExonGrouper
  
  attr_accessor :path_to_file
  attr_accessor :path_to_alignment
  attr_accessor :percent
  attr_accessor :organism_number

  attr_accessor :blossum_matrix

  attr_accessor :organisms
  attr_accessor :max_group
  attr_accessor :output_filename
  attr_accessor :exons_hash
  attr_accessor :local_borders

  def initialize(options = {})

    self.path_to_file = options["file"]
    self.path_to_alignment = options["alignment"]
    self.percent = options["percent"].to_i
    self.organism_number = options["organism_number"].to_i
    self.output_filename = options["output_filename"]

    self.blossum_matrix = DataProcessor.parse_blossum
    self.organisms = []
    self.exons_hash = {}
    self.local_borders = Hash.new { |hash, key| hash[key] = Array.new() }
  end

  def prepare_data
    # отбираются только первые self.organism_number организмов
    all_organisms = DataProcessor.new(self.path_to_file, self.path_to_alignment).prepare
    #self.organisms = all_organisms[0..(self.organism_number-1)]
    self.organisms = [all_organisms[0], all_organisms[17]]
    clear_output_file
    Organism.set_headers(output_filename)
    self.organisms.each{ |org| org.save_references(output_filename) }
  end

  def group
    time1 = Time.now
    count_cat_cat_statistics
  end

  def count_cat_cat_statistics
    # подготовка файла статистики
    CatCatMatcher.csv_header(output_filename)
    CatCatMatcher.clear_output_file(output_filename)
    pair_counter = 0
    
    self.organisms.each_with_index do |organism, index|
      puts "ORG: #{organism.name}"
      break if index+1 == self.organisms.length
      self.organisms[(index+1)..-1].each do |match_organism|
        splits = get_cat_cat_splits(organism, match_organism) # {"coords", "organism", "match_organism"}
        puts "  -> #{match_organism.name}"
        current_split_offset = 0
        exons_in_splits = ExonsSplitsIds.new(splits, organism, match_organism)
        splits["organism"].each_with_index do |org_part, part_index|

          coords = part_index == 0 ? [0, splits["coords"][0]] : [splits["coords"][part_index-1]+2, splits["coords"][part_index]] # +2, because of UU
          
          if part_index == (splits["organism"].length-1)
            coords = [splits["coords"][part_index-1], splits["coords"][part_index-1] + splits["organism"][-1].length + 2 ]
          end
          match_org_part = splits["match_organism"][part_index]
          
          sequences = [org_part, match_org_part]
          cat_cat_proxy = CatCatProxy.new(sequences, coords, self.blossum_matrix, organism, match_organism, 
                                          get_pair_id( part_index, organism.number, match_organism.number ),
                                          exons_in_splits, part_index, output_filename)
          cat_cat_matcher = CatCatMatcher.new(cat_cat_proxy)
          cat_cat_matcher.count_statistics
          cat_cat_matcher.save_to_csv(output_filename)
          # puts "O : #{org_part}"
          # puts "M : #{match_org_part}"
          # puts "\n"
          seq_1_bord = cat_cat_matcher.cat_cat_result.local_borders_seq_1
          seq_2_bord = cat_cat_matcher.cat_cat_result.local_borders_seq_2
          
          seq_1_bord = seq_1_bord.each_slice(2).to_a.map { |slice| slice.map{ |val| val+current_split_offset } }
          seq_2_bord = seq_2_bord.each_slice(2).to_a.map { |slice| slice.map{ |val| val+current_split_offset } }
          # puts "LB: #{seq_1_bord}"
          # puts "LB: #{seq_2_bord}"
          # res_1 = "C : "
          # seq_1_bord.each do |slice|
          #   res_1 += "#{organism.alignment[slice[0]..(slice[1])]}"
          # end
          # puts res_1

          # res_2 = "C : "
          # seq_2_bord.each do |slice|
          #   res_2 += "#{match_organism.alignment[slice[0]..(slice[1])]}"
          # end
          # puts res_2
          # puts "---------------"
          # if cat_cat_matcher.has_multiple?
          #   con_cat_matcher = ConCatMatcher.new(cat_cat_proxy)
          #   con_cat_matcher.align
          # end
          self.local_borders[organism.name] += seq_1_bord
          self.local_borders[match_organism.name] += seq_2_bord
          pair_counter += 1
          current_split_offset += org_part.length+2
        end
      end
      # puts "#{local_borders}"
      border_data = "organism:   #{organism.name}\n"
      border_data +="starts:     #{local_borders[organism.name].map{ |pair| pair[0] }}\n"
      border_data +="ends:       #{local_borders[organism.name].map{ |pair| pair[1] }}\n"

      File.open("#{self.output_filename}_borders.txt", 'a') { |file| file.write(border_data) }
    end
  end

  def get_cat_cat_splits(organism, match_organism)
    organism_coords = get_uu_coords(organism)
    match_organism_coords = get_uu_coords(match_organism)
    match_coords = organism_coords & match_organism_coords

    organism_parts = []
    match_organism_parts = []
    match_coords.each_with_index do |coord, index|
      if index == 0
        organism_parts << organism.alignment[0..(coord-1)]#.gsub("UU","--")
        match_organism_parts << match_organism.alignment[0..(coord-1)]#.gsub("UU","--")
        next
      end
      organism_parts << organism.alignment[(match_coords[index-1]+2)..(coord-1)]#.gsub("UU","--")
      match_organism_parts << match_organism.alignment[(match_coords[index-1]+2)..(coord-1)]#.gsub("UU","--")
    end
    if match_coords.length > 0
      organism_parts << organism.alignment[(match_coords[-1]+2)..(-1)]#.gsub("UU","--")
      match_organism_parts << match_organism.alignment[(match_coords[-1]+2)..(-1)]#.gsub("UU","--")
    end
    return ({"coords" => match_coords, "organism" => organism_parts, "match_organism" => match_organism_parts})
  end

  def get_exons_numbers_from_splits(splits)
    exons_organism_ids = []
  end

  def get_uu_coords(organism)
    i = -1
    all = []
    while i = organism.alignment.index('UU', i+1)
      all << i
    end
    return all
  end

  def draw_as_svg_rectangels(data_to_show)
    svg_width = self.organisms.first.alignment_length*2 + 200
    svg_height = self.organisms.count * 40 + 40
    output_file_name = self.output_filename + "_#{data_to_show}"
    File.open("#{output_file_name}.svg", 'w') do |file|
      file.write("<svg width=\"#{svg_width+100}\" height=\"#{svg_height}\">")
      file.write("<rect x=\"0\" y=\"0\" width=\"#{svg_width+100}\" height=\"#{svg_height}\" style=\"fill:white;\" />")
      self.organisms.each_with_index do |organism, index|
        draw_organism_line(index, svg_width, file)
        print_organism_name(organism, index, file)
        organism.exons.each do |exon|
          draw_exon_box(organism, index, exon, file, exon.send(data_to_show))
        end
      end
      #draw_exon_limits(file, svg_height)
      file.write("</svg>")
    end

    puts "exon number : #{organisms.map{ |org| org.exons.length }.inject(:+)}"
    `inkscape -z -e #{output_file_name}.png -w #{svg_width} -h #{svg_height} #{output_file_name}.svg`
  end
  
  def get_borders(aligner, exon, match_exon, sequence_data)
    data1 = [exon.uuid, match_exon.uuid, sequence_data[:pair_id],aligner.local_data['start_position_1'],aligner.local_data['end_position_1'],aligner.local_data['start_position_2'],aligner.local_data['end_position_2']]
    data1 = data1.join(",") + "\n"
    data2 = [match_exon.uuid, exon.uuid, sequence_data[:pair_id],aligner.local_data['start_position_2'],aligner.local_data['end_position_2'],aligner.local_data['start_position_1'],aligner.local_data['end_position_1']]
    data2 = data2.join(",") + "\n"
    return data1+data2
  end

  def draw_exon_limits(file, svg_height)
    starts = Set.new
    finish = Set.new
    exons = self.organisms.map(&:exons).flatten
    exons.each do |exon|
      starts.add(exon.start)
      finish.add(exon.finish)
    end
    starts.each { |st| draw_border_line(st, "rgb(232, 18, 18)", file, svg_height) }
    finish.each { |fn| draw_border_line(fn, "rgb(18, 18, 232)", file, svg_height) }
  end

private

  def clear_output_file
    header = "ex1,ex2,pair_id,st1,end1,st2,end2\n"
    File.open("#{self.output_filename}_borders.txt", 'w') { |file| file.write(header) }
    IterSaver.header(self.output_filename)
    #ConCatMatcher.header(self.output_filename)
    #ConCatMatcher.alignements_header(self.output_filename)
  end

  def get_pair_id(number, org_num, match_org_num)
    return "A"+"0"*(3-number.to_s.length) + number.to_s + "_" + 
           "0"*(3-org_num.to_s.length) + org_num.to_s + "_" + 
           "0"*(3-match_org_num.to_s.length) + match_org_num.to_s
  end

  def print_organism_name(organism, index, file)
    y_coord = 40*(index+1) - 10
    draw_text(10, y_coord, organism.name, file)
  end

  def draw_organism_line(index, svg_width, file)
    y_coord = 40*(index+1)
    x_start_coords = 100
    x_end_coords = svg_width - 100
    draw_line(x_start_coords, y_coord, x_end_coords, y_coord, "rgb(0,0,0)", file)
  end

  def draw_border_line(x_start, color, file, svg_height)
    x1 = (x_start + 100)*2
    y1 = 0
    draw_line(x1, y1, x1, svg_height, color, file)
  end

  def draw_exon_box(organism, index, exon, file, data_to_show)
    y_coord = 40*(index+1)
    x_start_coords = 100
    width = exon.finish - exon.start
    color = exon.get_svg_color

    x = (exon.start+x_start_coords)*2
    y = y_coord-15
    width = width*2
    height = 30
    draw_rect(x, y, width, height, color, file)

    x_text = (exon.start+x_start_coords)*2
    draw_text(x_text, y_coord, data_to_show, file)
  end

  def draw_line(x1,y1,x2,y2,color,file)
    file.write("<line x1=\"#{x1}\" y1=\"#{y1}\" x2=\"#{x2}\" y2=\"#{y2}\" style=\"stroke:#{color};stroke-width:1\" />")
  end

  def draw_rect(x,y,width,height,color,file)
    file.write("<rect x=\"#{x}\" y=\"#{y}\" width=\"#{width}\" height=\"#{height}\" style=\"fill:#{color};stroke:black;stroke-width:0\" />\n")
  end

  def draw_text(x,y,text,file)
    file.write("<text x=\"#{x}\" y=\"#{y}\" fill=\"black\">#{text}</text>")
  end

end