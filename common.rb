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
      @ids_for_organism = count_ids(splits["organism"], organism_number)
    end
    return @ids_for_organism 
  end

  def get_ids_for_match_organism
    if @ids_for_match_organism.nil?
      @ids_for_match_organism = count_ids(splits["match_organism"], match_organism_number)
    end
    return @ids_for_match_organism
  end

  def count_ids(org_splits, org_number)
    exon_org_ids = []
    last_id_count = 0
    org_splits.each do |org_split|
      number_of_exons_in_split = org_split.split("UU").length
      current_ids = (1..( number_of_exons_in_split )).to_a.map{ |num|  "#{ (org_number+1)*100+num+last_id_count}"  }
      last_id_count += number_of_exons_in_split
      exon_org_ids << current_ids
    end
    return exon_org_ids
  end

end

class CatCatProxy

  attr_accessor :sequences
  attr_accessor :coords
  attr_accessor :blosum
  attr_accessor :organism
  attr_accessor :match_organism
  attr_accessor :pair_id
  attr_accessor :exons_ids
  attr_accessor :split_index

  def initialize(sequences = [], coords = [], blosum = nil, organism, match_organism, pair_id, exons_ids, split_index)
    self.sequences = sequences
    self.coords = coords
    self.blosum = blosum.nil? ? DataProcessor.parse_blossum : blosum
    self.organism = organism
    self.match_organism = match_organism
    self.pair_id = pair_id
    self.exons_ids = exons_ids
    self.split_index = split_index
  end

  def seq_1
    self.sequences[0]
  end

  def seq_2 
    self.sequences[1]
  end

  def get_org_exons 
    exons_ids.get_ids_for_organism[self.split_index].join(",")
  end

  def get_match_org_exons
    exons_ids.get_ids_for_match_organism[self.split_index].join(",")
  end

  def get_org_exons_count 
    exons_ids.get_ids_for_organism[self.split_index].length
  end

  def get_match_org_exons_count
    exons_ids.get_ids_for_match_organism[self.split_index].length
  end

end

class CatCatResult

  attr_accessor :added_spaces
  attr_accessor :matching_spaces
  attr_accessor :matching_letters
  attr_accessor :mismatching_letters
  attr_accessor :affine_score
  attr_accessor :usual_score

  attr_accessor :deletions_1
  attr_accessor :deletions_2

  attr_accessor :seq_1_score
  attr_accessor :seq_2_score

  attr_accessor :local_score
  attr_accessor :local_seq_1
  attr_accessor :local_seq_2
  attr_accessor :local_self_1_score
  attr_accessor :local_self_2_score

  attr_accessor :local_borders_seq_1
  attr_accessor :local_borders_seq_2

  def initialize
    self.added_spaces = 0
    self.matching_spaces = 0
    self.matching_letters = 0
    self.mismatching_letters = 0
    self.affine_score = 0
    self.usual_score = 0

    self.deletions_1 = 0
    self.deletions_2 = 0

    self.seq_1_score = 0
    self.seq_2_score = 0

    self.local_score = 0
    self.local_seq_1 = ""
    self.local_seq_2 = ""
    self.local_self_1_score = 0
    self.local_self_2_score = 0

    self.local_borders_seq_1 = []
    self.local_borders_seq_2 = []
  end

  def check
    if ((seq_1_score < usual_score) || (seq_2_score < usual_score))
      raise "self score > usual score"  
    end

  end

  def print
    puts ": #{added_spaces}"
    puts ": #{matching_spaces}"
    puts ": #{matching_letters}"
    puts ": #{mismatching_letters}"
    puts ": #{affine_score}"
    puts ": #{usual_score}"
    puts ": #{deletions_1}"
    puts ": #{deletions_2}"
    puts ": #{seq_1_score}"
    puts ": #{seq_2_score}"
  end

  def save(file_path, proxy)
    aff_seq_1 = (self.affine_score/self.seq_1_score.to_f).round(2)
    aff_seq_2 = (self.affine_score/self.seq_2_score.to_f).round(2)
    rloc_1 = (self.local_score+0.01)/aff_seq_1.to_f
    rloc_2 = (self.local_score+0.01)/aff_seq_2.to_f
    data = [
      proxy.pair_id,
      proxy.organism.name,
      proxy.match_organism.name,
      proxy.organism.code_name,
      proxy.match_organism.code_name,
      proxy.get_org_exons,
      proxy.get_match_org_exons,
      proxy.get_org_exons_count,
      proxy.get_match_org_exons_count,
      [proxy.get_org_exons_count,proxy.get_match_org_exons_count].max,
      [proxy.get_org_exons_count,proxy.get_match_org_exons_count].min,
      proxy.coords[0],
      proxy.coords[1],
      proxy.coords[1]-proxy.coords[0],
      self.seq_1_score,
      self.seq_2_score,
      self.usual_score,
      self.affine_score,
      aff_seq_1,
      aff_seq_2,
      [aff_seq_1, aff_seq_2].min,
      "",
      ""
    ]
    #puts self.affine_score
    #local values
    data += [
      self.local_score,
      self.local_score.to_f-self.affine_score,
      self.local_self_1_score,
      self.local_self_2_score,
      self.local_score.to_f/[self.local_self_1_score,self.local_self_2_score].max,
      self.local_score/self.local_self_1_score.to_f,
      self.local_score/self.local_self_2_score.to_f,
      "",
      "",
      [self.local_self_1_score, self.local_self_2_score].min.to_f,
      rloc_1,
      rloc_2,
      [rloc_1,rloc_2].min,
      [rloc_1,rloc_2].max,
    ]
    CSV.open("#{file_path}.csv", "a") { |csv| csv << data}
    save_allignement(file_path, proxy)
  end

  def save_allignement(file_path, proxy)
    output = ""
    output += "#{proxy.pair_id} Coords: #{proxy.coords} Exons_1: #{proxy.get_org_exons} Exons_2: #{proxy.get_match_org_exons}\n"
    output += "\"#{proxy.seq_1}\"\n\"#{proxy.seq_2}\"\n"
    File.open("#{file_path}_allignements.txt", 'a') { |file| file.write(output) }

    output = ""
    output += "#{proxy.pair_id} Coords: #{proxy.coords} Exons_1: #{proxy.get_org_exons} Exons_2: #{proxy.get_match_org_exons}\n"
    output += "\"#{self.local_seq_1}\"\n\"#{self.local_seq_2}\"\n"
    File.open("#{file_path}_local_allignements.txt", 'a') { |file| file.write(output) }
  end

  def save_references(file_path, proxy)

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