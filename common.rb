class Split

  SINGLE = 1
  NORMAL = 2

  attr_accessor :organism
  attr_accessor :match_organism
  attr_accessor :seq_1
  attr_accessor :seq_2

  attr_accessor :start_coord
  attr_accessor :end_coord

  attr_accessor :type
  attr_accessor :org_exons_1
  attr_accessor :org_exons_2

  def initialize(organism, match_organism, seq_1="", seq_2="", start_coord, end_coord, type)
    self.organism = organism
    self.match_organism = match_organism
    self.seq_1 = seq_1
    self.seq_2 = seq_2
    self.start_coord = start_coord
    self.end_coord = end_coord
    self.type = type
    self.org_exons_1 = []
    self.org_exons_2 = []
  end

  def print
    puts "S: #{start_coord}"
    puts "E: #{end_coord}"
    puts "E1: #{org_exons_1}"
    puts "E2: #{org_exons_2}"
    puts "#{seq_1}"
    puts "#{seq_2}"
  end

end

class ExonsSplitsIds

  attr_accessor :organism_number
  attr_accessor :match_organism_number
  attr_accessor :splits

  def initialize(splits, organism, match_organism)
    self.organism_number = organism.number
    self.match_organism_number = match_organism.number
    self.splits = splits
  end

  def get_ids_for_organism
    if @ids_for_organism.nil?
      @ids_for_organism = count_ids(splits.map(&:seq_1), organism_number)
    end
    return @ids_for_organism 
  end

  def get_ids_for_match_organism
    if @ids_for_match_organism.nil?
      @ids_for_match_organism = count_ids(splits.map(&:seq_2), match_organism_number)
    end
    return @ids_for_match_organism
  end

  def count_ids(org_splits, org_number)
    exon_org_ids = []
    last_id_count = 0
    org_splits.each do |org_split|
      number_of_exons_in_split = org_split.gsub("-","").split("UU").length
      current_ids = (1..( number_of_exons_in_split )).to_a.map{ |num|  "#{ (org_number+1)*100+num+last_id_count}"  }
      last_id_count += number_of_exons_in_split
      exon_org_ids << current_ids
    end
    return exon_org_ids
  end

  def set_to_splits
    splits.each_with_index do |split, index|
    split.org_exons_1 = get_ids_for_organism[index]
    split.org_exons_2 = get_ids_for_match_organism[index]
    end
  end

end

class CatCatProxy

  attr_accessor :pair_id
  attr_accessor :split
  attr_accessor :blosum

  def initialize(pair_id, split, blosum )
    self.pair_id = pair_id
    self.split = split
    self.blosum = blosum
  end

end

class CatCatResult

  attr_accessor :local_iters
  attr_accessor :locals

  def initialize
    self.local_iters = 0
  end

  def check
    if ((seq_1_score < usual_score) || (seq_2_score < usual_score))
      raise "self score > usual score"  
    end

  end


  def local_seq_1
    return locals.map{ |local| local.get_formatted_seq(1) }.join("")
  end

  def local_seq_2
    return locals.map{ |local| local.get_formatted_seq(2) }.join("")
  end

  def local_seq_1_for_coords
    return locals.map{ |local| local.get_for_coords_seq(1) }.join("")
  end

  def local_seq_2_for_coords
    return locals.map{ |local| local.get_for_coords_seq(2) }.join("")
  end

  def save_alignment(file_path, proxy)
    output = ""
    output += "#{proxy.pair_id} Coords: #{proxy.coords} Exons_1: #{proxy.get_org_exons} Exons_2: #{proxy.get_match_org_exons}\n"
    output += "\"#{proxy.seq_1}\"\n\"#{proxy.seq_2}\"\n"
    File.open("#{file_path}_alignments.txt", 'a') { |file| file.write(output) }

    output = ""
    output += "#{proxy.pair_id} Coords: #{proxy.coords} Exons_1: #{proxy.get_org_exons} Exons_2: #{proxy.get_match_org_exons}\n"
    output += "\"#{self.local_seq_1}\"\n\"#{self.local_seq_2}\"\n"
    File.open("#{file_path}_local_alignments.txt", 'a') { |file| file.write(output) }
  end

end

class Points

  attr_accessor :organism
  attr_accessor :match_organism
  attr_accessor :points

  def initialize(organism, match_organism)
    self.organism = organism
    self.match_organism = match_organism
    self.points = []
  end

  def make_points(cat_res, split)
    local_seq_1 = cat_res.local_seq_1_for_coords
    local_seq_2 = cat_res.local_seq_2_for_coords
    good_uu_coords_seq_1 = get_good_uu_coords(local_seq_1)
    good_uu_coords_seq_2 = get_good_uu_coords(local_seq_2)
    bad_uu_coords_seq_1 = get_bad_uu_coords(local_seq_1)
    bad_uu_coords_seq_2 = get_bad_uu_coords(local_seq_2)
    
    puts "#{split.seq_1}"
    puts "#{split.seq_2}"
    puts "----"
    puts "#{local_seq_1}"
    puts "#{local_seq_2}"
    puts "#{good_uu_coords_seq_1}"
    puts "#{good_uu_coords_seq_2}"
    puts "#{bad_uu_coords_seq_1}"
    puts "#{bad_uu_coords_seq_2}"
    puts "--------------------------"

    good_uu_coords_seq_1.each do |coord|
      puts "Y : #{split.seq_2[0..(get_coords_in_align(split.seq_2, local_seq_2, coord))]}"
    end
    bad_uu_coords_seq_1.each do |coord|
      puts "Y : #{split.seq_2[0..(get_coords_in_align(split.seq_2, local_seq_2, coord))]} "
    end

    good_uu_coords_seq_2.each do |coord|
      puts "Y : #{split.seq_1[0..(get_coords_in_align(split.seq_1, local_seq_1, coord))]} "
    end
    bad_uu_coords_seq_2.each do |coord|
      puts "Y : #{split.seq_1[0..(get_coords_in_align(split.seq_1, local_seq_1, coord))]} "
    end



  end

  private

  def get_good_uu_coords(sequence)
    i = -1
    all = []
    while i = sequence.index('UU', i+1)
      all << i
    end
    return all
  end

  def get_bad_uu_coords(sequence)
    i = -1
    all = []
    while i = sequence.index('uu', i+1)
      all << i
    end
    return all
  end

  def get_coords_in_align(sequence, local, coord)
    return (get_align_coords(sequence, get_letters_before(local,coord))-1)
  end

  def get_letters_before(sequence, coord)
    return sequence[0..(coord-1)].gsub("-","").length
  end

  def get_align_coords(sequence, coord)
    result = 0
    (sequence.length-1).times do |pos|
      coord -= 1 if sequence[pos] != "-"
      #result << sequence[pos]
      result += 1
      break if coord == 0
    end
    return result
  end


end

class SplitPoint

  OCCURATE = 0
  FUZZY = 1
  UU = 2

  attr_accessor :start
  attr_accessor :finish
  attr_accessor :type
  attr_accessor :seq_number

  def close_reg(finish)
    if type == SplitPoint::FUZZY
      self.finish = finish
    else
      self.finish = start
    end
  end

  def is_occurate?
    return type == SplitPoint::OCCURATE
  end

  def is_fuzzy?
    return type == SplitPoint::FUZZY
  end

  def is_uu?
    return type == SplitPoint::UU
  end

end

class LocalResult

  attr_accessor :start_positions
  attr_accessor :end_positions
  attr_accessor :align_1
  attr_accessor :align_2
  attr_accessor :score

  def initialize(results)
    self.start_positions = results["start_positions"]
    self.end_positions = results["end_positions"]
    self.align_1 = results["align_1"]
    self.align_2 = results["align_2"]
    self.score = results["score"]
    #File.open("PAIRSTAT_fras_30(ALL)", "a") { |file| file.write("#{self.score},") }
  end

  def from_start_1?
    return start_1 == 0
  end

  def from_start_2?
    return start_2 == 0
  end

  def till_end_1?(seq)
    return end_1 == (seq.length-1)
  end

  def till_end_2?(seq)
    return end_2 == (seq.length-1)    
  end

  def start_1
    return self.start_positions[0]
  end

  def start_2
    return self.start_positions[1]
  end

  def end_1
    return self.end_positions[0]
  end

  def end_2
    return self.end_positions[1]
  end

  def aligns
    return [align_1, align_2]
  end
end

class LocalReqursiveResult

  GOOD = 1
  BAD = 0

  attr_accessor :seq_1
  attr_accessor :seq_2
  attr_accessor :index
  attr_accessor :type
  attr_accessor :coords_org_1
  attr_accessor :coords_org_2

  def initialize(seq_1, seq_2, type, index = 9999)
    self.seq_1 = seq_1
    self.seq_2 = seq_2
    self.index = index
    self.type = type
  end

  def get_formatted
    return [get_formatted_seq(1), get_formatted_seq(2)]
  end

  def get_for_coords
    return [get_for_coords_seq(1), get_for_coords_seq(2)]
  end

  def get_formatted_seq(seq_number)
    if type == GOOD
      return "+#{index}+"+self.send("seq_#{seq_number}")+"+#{index}+"
    else
      return "|"+self.send("seq_#{seq_number}")+"|"
    end
  end

  def get_for_coords_seq(seq_number)
    return self.send("seq_#{seq_number}") if type == GOOD
    return self.send("seq_#{seq_number}").downcase
  end

  def get_raw
    return [seq_1, seq_2]
  end
end

class IterSaver

  def self.save(proxy, locals)
    iters_number = locals.map(&:index).delete_if { |val| val == 9999 }.max
    locals.sort_by(&:index).each do |local|
      next if local.type == LocalReqursiveResult::BAD
      result = []
      result << proxy.pair_id
      result << proxy.organism.code_name
      result << proxy.match_organism.code_name
      result << local.seq_1.length
      result << iters_number
      result << local.index
      result << local.seq1_exons_ids.join(",")
      result << local.seq2_exons_ids.join(",")
      result << local.seq1_exons_ids.length
      result << local.seq2_exons_ids.length
      result << local.score
      result << local.score_seq_1
      result << local.score_seq_2
      result << (local.score.to_f/local.score_seq_1).round(2)
      result << (local.score.to_f/local.score_seq_2).round(2)
      CSV.open("#{proxy.file_path}_iters.csv", "a") { |csv| csv << result } 
    end

  end

  def self.header(file_path)
    CSV.open("#{file_path}_iters.csv", "w") do |csv|
      header = ["pair_id",
                "spec_1",
                "spec_2",
                "length",
                "num_iter",
                "iter",
                "exons_ids_1",
                "exons_ids_2",
                "size_ex_1",
                "size_ex_2",
                "local_score",
                "self_local_score_1",
                "self_local_score_2",
                "local_score_1_coef",
                "local_score_2_coef",
            ]
      csv << header
    end
  end

end