require_relative "common.rb"
require_relative "organism.rb"
require_relative "cat_cat_local_aligner.rb"

class ConCatMatcher

  attr_accessor :proxy

  def initialize(proxy)
    self.proxy = proxy
  end


  def align
    no_gaps_seqs_1 = proxy.seq_1.gsub("-","")
    no_gaps_seqs_2 = proxy.seq_2.gsub("-","")

    exons_coords_seq_1 = get_exons_coords(proxy.seq_1)
    exons_coords_seq_2 = get_exons_coords(proxy.seq_2)

    locals_seq_1 = []
    locals_seq_2 = []

    no_gaps_seqs_1.split("UU").each do |seq_part|
      locals_seq_1 << CatCatLocalAligner.new(seq_part, proxy.seq_2, proxy.blosum).align
    end

    no_gaps_seqs_2.split("UU").each do |seq_part|
      locals_seq_2 << CatCatLocalAligner.new(seq_part, proxy.seq_1, proxy.blosum).align
    end
      
    # puts "#{proxy.seq_1}"
    # puts "C #{exons_coords_seq_1}"

    # locals_seq_1.each do |local|
    #   puts "#{local.align_1}"
    #   puts "#{local.align_2}"
    #   puts "#{local.start_positions[1]} : #{local.end_positions[1]}"
    #   puts "--"
    # end

    # puts "#{proxy.seq_2}"
    # puts "C #{exons_coords_seq_2}"
    # locals_seq_2.each do |local|
    #   puts "#{local.align_1}"
    #   puts "#{local.align_2}"
    #   puts "#{local.start_positions[1]} : #{local.end_positions[1]}"
    #   puts "--"
    # end


    locals_seq_1.each_with_index do |local, local_index|
      local.current_exon = self.proxy.get_org_exons_raw[local_index]
      connected_exons_indexes = get_connected_exons(local, exons_coords_seq_2) 
      local.connected_exons = connected_exons_indexes.map { |val| val = self.proxy.get_match_org_exons_raw[val] }
      local.connected_coords = connected_exons_indexes.map { |val| exons_coords_seq_2[val] }
      local.first_coverage, local.last_coverage = get_coverage(local, local.connected_coords)
    end
    locals_seq_2.each_with_index do |local, local_index|
      local.current_exon = self.proxy.get_match_org_exons_raw[local_index]
      connected_exons_indexes = get_connected_exons(local, exons_coords_seq_1)
      local.connected_exons = connected_exons_indexes.map { |val| val = self.proxy.get_org_exons_raw[val] }
      local.connected_coords = connected_exons_indexes.map { |val| exons_coords_seq_1[val] }
      local.first_coverage, local.last_coverage = get_coverage(local, local.connected_coords)
    end

    locals_seq_1.each { |local| save(local) }
    locals_seq_2.each { |local| save(local) }

  end

  def get_coverage(local, coords)
    first_coverage, last_coverage = 0, 0 
    if local.end_positions[1] > coords[0][1]
      first_coverage = (coords[0][1] - local.start_positions[1]).to_f/(coords[0][1] - coords[0][0])
    end
    if local.end_positions[1] <= coords[0][1]
      first_coverage = (local.end_positions[1] - local.start_positions[1]).to_f/(coords[0][1] - coords[0][0])
    end
    if local.start_positions[1] < coords[-1][0]
      last_coverage = (local.end_positions[1] - coords[-1][0]).to_f/(coords[-1][1] - coords[-1][0])
    end
    if local.start_positions[1] >= coords[-1][0]
      last_coverage = (local.end_positions[1] - local.start_positions[1]).to_f/(coords[-1][1] - coords[-1][0])
    end
    return [first_coverage, last_coverage]
  end

  def get_connected_exons(local, exons_coords)
    exons = []
    exons_coords.each_with_index do |coord, index|
      if (local.start_positions[1] >= coord[0]) & (local.start_positions[1] <= coord[1])
        exons << index
      end
      if (local.end_positions[1] >= coord[0]) & (local.end_positions[1] <= coord[1])
        exons << index
      end
    end
    return exons.uniq
  end
  
  def get_exons_coords(seq)
    coords = []
    uu_coords = get_uu_coords(seq)
    if uu_coords.length > 0
      coords << [0,uu_coords[0]-1]
      uu_coords.each_with_index do |coord, index|
        break if index == uu_coords.length-1
        coords << [coord+2, uu_coords[index+1]-1]
      end
      coords << [uu_coords[-1]+2, seq.length-1]
    else
      coords = [[0,seq.length-1]]
    end 
    return coords
  end

  def get_uu_coords(seq)
    i = -1
    all = []
    while i = seq.index('UU', i+1)
      all << i
    end
    return all
  end

  def save(local)
    scorer = CatCatLocalAligner.new("", "", proxy.blosum)
    self_1_score = scorer.count_score(local.align_1, local.align_1).to_f
    self_2_score = scorer.count_score(local.align_2, local.align_2).to_f
    data = [
            proxy.pair_id,
            proxy.organism.code_name,
            proxy.match_organism.code_name,
            proxy.get_org_exons,
            proxy.get_match_org_exons,
            proxy.get_org_exons_raw.length,
            proxy.get_match_org_exons_raw.length,
            local.current_exon,
            local.connected_exons[0],
            local.connected_exons[-1],
            local.first_coverage.round(2),
            local.last_coverage.round(2),
            local.score,
            self_1_score,
            self_2_score,
            (local.score/self_1_score).round(2),
            (local.score/self_2_score).round(2)
    ]
    CSV.open("#{proxy.file_path}_con_cat_statistics.csv", 'a') { |csv| csv << data}
    alignement_data = ""
    alignement_data += "#{proxy.pair_id} Connected exons coords: #{local.connected_coords} Exon: #{local.current_exon} Connected: #{local.connected_exons.join(',')} \n"
    alignement_data += "#{local.align_1}\n#{local.align_2}\n"
    File.open("#{proxy.file_path}_con_cat_statistics_alignements.txt", 'a') { |file| file.write(alignement_data) }
  end

  def self.alignements_header(output_filename)
    File.open("#{output_filename}_con_cat_statistics_alignements.txt", 'w') { |file| file.write("") }
  end

  def self.header(output_filename)
    header = [
              "pair_id",
              "Spec1",
              "Spec2",
              "exons_org_1",
              "exons_org_2",
              "#exons_org_1",
              "#exons_org_2",
              "current_exon",
              "first_domain_index",
              "last_domain_index",
              "first_cover_coef",
              "last_cover_coef",
              "local_score",
              "self_1_score",
              "self_2_score",
              "self_1_coef",
              "self_2_coef",
              "\n"
    ]
    File.open("#{output_filename}_con_cat_statistics.csv", 'w') { |file| file.write(header.join(",")) }
  end




end