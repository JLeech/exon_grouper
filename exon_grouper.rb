require "csv"
require "json"
require "set"

require_relative "exon.rb"
require_relative "organism.rb"
require_relative "data_processor"
require_relative "cat_cat_matcher.rb"
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
    self.organisms = all_organisms[0..(self.organism_number-1)]
    #self.organisms = [all_organisms[0],all_organisms[2]]
    #self.organisms = all_organisms[0..10]
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
    #CatCatMatcher.csv_header(output_filename)
    CatCatMatcher.clear_output_file(output_filename)
    pair_counter = 0
    file = File.open("/home/eve/Documents/biology/exon_grouper/results/result_dscam_uuid.svg", 'a')
    self.organisms.each_with_index do |organism, index|
      next if organism.name == "Colius_striatus"
      puts "ORG: #{organism.name}"
      puts "EXN: #{organism.exons.map(&:get_coords)}"
      match_org_index = index
      break if index+1 == self.organisms.length
      self.organisms[(index+1)..-1].each do |match_organism|
        match_org_index += 1
        splits = get_cat_cat_splits(organism, match_organism) # {"coords", "organism", "match_organism"}
        exons_in_splits = ExonsSplitsIds.new(splits, organism, match_organism)
        exons_in_splits.set_to_splits
        splits.each_with_index do |split, split_index|
          split_id = get_split_id( split_index, organism.number, match_organism.number )
          split.split_id = split_id
          proxy = CatCatProxy.new(split_id, split, self.blossum_matrix)
          cat_cat_matcher = CatCatMatcher.new(proxy)
          cat_cat_result = cat_cat_matcher.make_local
          
          split.points = Points.new(organism, match_organism)
          split.points.make_points(cat_cat_result, split)

          split.for_split_table(output_filename)
          cat_cat_result.save_alignment(output_filename, proxy)

          split.for_point_table(output_filename)
          split.points.organism_uniq_points.each do |point|
              draw_box(index, point,file, "")
          end
          split.points.match_organism_uniq_points.each do |point|
              draw_box(match_org_index, point,file, "")
          end
        end
      end
    end
    file.write("</svg>")
  end

  def get_cat_cat_splits(organism, match_organism)
    organism_coords = get_uu_coords(organism)
    match_organism_coords = get_uu_coords(match_organism)
    aligned_coords = organism_coords & match_organism_coords

    splits = []

    aligned_coords.each_with_index do |coord, index|
      if index == 0
        organism_part = organism.alignment[0..(coord-1)]#.gsub("UU","--")
        match_organism_part = match_organism.alignment[0..(coord-1)]#.gsub("UU","--")
        start_coord = 0
        end_coord = coord-1
        splits << Split.new(organism, match_organism, organism_part, match_organism_part, start_coord,end_coord,Split::SINGLE)
        next
      end
      organism_part = organism.alignment[(aligned_coords[index-1]+2)..(coord-1)]#.gsub("UU","--")
      match_organism_part = match_organism.alignment[(aligned_coords[index-1]+2)..(coord-1)]#.gsub("UU","--")
      start_coord = aligned_coords[index-1]+2
      end_coord = coord-1
      splits << Split.new(organism, match_organism, organism_part, match_organism_part, start_coord,end_coord,Split::NORMAL)
    end

    organism_part = organism.alignment[(aligned_coords[-1]+2)..(-1)]#.gsub("UU","--")
    match_organism_part = match_organism.alignment[(aligned_coords[-1]+2)..(-1)]#.gsub("UU","--")      
    start_coord = aligned_coords[-1]+2
    end_coord = organism.alignment.length-1
    splits << Split.new(organism, match_organism, organism_part, match_organism_part, start_coord,end_coord,Split::SINGLE)

    return splits
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
          draw_box(index, exon, file, exon.send(data_to_show))
        end
      end
      #draw_exon_limits(file, svg_height)
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
    # header = "ex1,ex2,pair_id,st1,end1,st2,end2\n"
    # File.open("#{self.output_filename}_borders.txt", 'w') { |file| file.write(header) }
    # IterSaver.header(self.output_filename)
    #Points.header(self.output_filename)
    Split.header(self.output_filename)
    Split.point_header(self.output_filename)
    #ConCatMatcher.header(self.output_filename)
    #ConCatMatcher.alignements_header(self.output_filename)
  end

  def get_split_id(number, org_num, match_org_num)
    return "A"+"0"*(3-number.to_s.length) + number.to_s + "_" + 
           "0"*(3-org_num.to_s.length) + org_num.to_s + "_" + 
           "0"*(3-match_org_num.to_s.length) + match_org_num.to_s
  end

  def uu_start_match(org_coords, match_org_coords, aligned_coords)
    return ((org_coords[0] == aligned_coords[0]) && (match_org_coords[0] == aligned_coords[0]))
  end

  def uu_end_match(org_coords, match_org_coords, aligned_coords)
    return ((org_coords[-1] == aligned_coords[-1]) && (match_org_coords[-1] == aligned_coords[-1]))
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

  def draw_box(index, box_data, file, data_to_show)
    if box_data.class == SplitPoint
      y_coord = 40*(index+1)-4
    else 
      y_coord = 40*(index+1)
    end
    x_start_coords = 100
    width = box_data.finish - box_data.start
    color = box_data.get_svg_color
    x = (box_data.start+x_start_coords)*2
    y = y_coord-15
    width = 1 if width == 0 

    width = width*2
    height = 30
    draw_rect(x, y, width, height, color, file)

    x_text = (box_data.start+x_start_coords)*2
    draw_text(x_text, y_coord, data_to_show, file)
  end

  def draw_line(x1,y1,x2,y2,color,file)
    file.write("<line x1=\"#{x1}\" y1=\"#{y1}\" x2=\"#{x2}\" y2=\"#{y2}\" style=\"stroke:#{color};stroke-width:1\" />")
  end

  def draw_rect(x,y,width,height,color,file)
    file.write("<rect x=\"#{x}\" y=\"#{y}\" width=\"#{width}\" height=\"#{height}\" style=\"fill:#{color};stroke:black;stroke-width:0;\" />\n")
  end

  def draw_text(x,y,text,file)
    file.write("<text x=\"#{x}\" y=\"#{y}\" fill=\"black\">#{text}</text>")
  end

end