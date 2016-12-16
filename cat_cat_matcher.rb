require_relative "exon.rb"
require_relative "data_processor.rb"
require_relative "cat_cat_local_aligner.rb"

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
      proxy.get_org_exons,
      proxy.get_match_org_exons,
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
      "",
      (self.affine_score.to_f/[seq_1_score, seq_2_score].max).round(2)
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
    output += "#{proxy.pair_id} #{proxy.coords} Exons_1: #{proxy.get_org_exons} Exons_2: #{proxy.get_match_org_exons}\n"
    output += "\"#{self.local_seq_1}\"\n\"#{self.local_seq_2}\"\n"
    File.open("#{file_path}_local_allignements.txt", 'a') { |file| file.write(output) }
  end


end

class CatCatMatcher

  attr_accessor :cat_cat_proxy
  attr_accessor :cat_cat_result

  def initialize(cat_cat_proxy)
    self.cat_cat_proxy = cat_cat_proxy
    self.cat_cat_result = CatCatResult.new
  end

  def count_statistics
    cat_cat_result.affine_score, cat_cat_result.usual_score = count_blosum(cat_cat_proxy.seq_1, cat_cat_proxy.seq_2, true)
    cat_cat_result.seq_1_score, _ = count_blosum(cat_cat_proxy.seq_1, cat_cat_proxy.seq_1)
    cat_cat_result.seq_2_score, _ = count_blosum(cat_cat_proxy.seq_2, cat_cat_proxy.seq_2)
    make_local
    cat_cat_result.check
  end

  def make_local
    no_gap_seq_1 = cat_cat_proxy.seq_1.gsub("-","")
    no_gap_seq_2 = cat_cat_proxy.seq_2.gsub("-","")
    locals = local_recursive(no_gap_seq_1, no_gap_seq_2)
    local_1 = ""
    local_2 = ""
    locals.each do |local_split|
      local_1 += local_split[0]
      local_2 += local_split[1]
    end
    cat_cat_result.local_seq_1 = local_1
    cat_cat_result.local_seq_2 = local_2
    cat_cat_result.local_score, _ = count_blosum(local_1, local_2)
    cat_cat_result.local_self_1_score, _ = count_blosum(local_1, local_1)
    cat_cat_result.local_self_2_score, _ = count_blosum(local_2, local_2)
    cat_cat_result.local_borders_seq_1 = get_local_borders(local_1, cat_cat_proxy.seq_1)
    cat_cat_result.local_borders_seq_2 = get_local_borders(local_2, cat_cat_proxy.seq_2)
    raise "local borders error" if cat_cat_result.local_borders_seq_1.length%2 != 0
    raise "local borders error" if cat_cat_result.local_borders_seq_2.length%2 != 0
    # puts "#{self.cat_cat_result.usual_score} : #{self.cat_cat_result.affine_score}"
    # puts "#{blos_loc}"
    # puts ""
    # puts "#{cat_cat_proxy.seq_1}"
    # puts "#{cat_cat_proxy.seq_2}"
    # puts "L1"
    # puts "#{local_1}"
    # puts "L2"
    # puts "#{local_2}"
    # puts "\n-----------\n"
  end

  def local_recursive(seq_1, seq_2, len_coef = 1.0)
    if stop_results?(seq_1, seq_2, len_coef)
      return [margin(seq_1,seq_2)]
    end
    local_result = CatCatLocalAligner.new(seq_1, seq_2, cat_cat_proxy.blosum).align
    min_seq_length = [seq_1, seq_2].map(&:length).min.to_f
    left_part_1 = ""
    left_part_2 = ""
    right_part_1 = seq_1[(local_result.end_1+1)..(-1)]
    right_part_2 = seq_2[(local_result.end_2+1)..(-1)]
    if local_result.start_1 > 0
      left_part_1 = seq_1[0..(local_result.start_1-1)]
    end
    if local_result.start_2 > 0
      left_part_2 = seq_2[0..(local_result.start_2-1)]
    end
    left_length_coef = [left_part_1,left_part_2].map(&:length).min.to_f/min_seq_length
    right_length_coef = [right_part_1,right_part_2].map(&:length).min.to_f/min_seq_length
    final = ""
    if (local_result.align_1.length < 5) || (local_result.align_2.length < 5)
      puts "final!"
      final = [margin(seq_1, seq_2)] 
    else
      final = local_recursive(left_part_1,left_part_2,left_length_coef) + [["+#{local_result.align_1}+", "+#{local_result.align_2}+"]] + local_recursive(right_part_1,right_part_2,right_length_coef)
    end
    return final
  end

  def stop_results?(seq_1, seq_2, len_coef)
    if (len_coef < 0.2) || (seq_1.length < 5) || (seq_2.length < 5)
      return true
    end
    return false
  end

  def margin(seq_1, seq_2)
    if seq_1.length > seq_2.length
      seq_2 += "-"*(seq_1.length-seq_2.length)
    elsif seq_1.length < seq_2.length
      seq_1 += "-"*(seq_2.length-seq_1.length)
    end
    return ["|"+seq_1+"|", "|"+seq_2+"|"]
  end

  def count_blosum(main_allignement, matching_allignement, global = false)
    current_affine_score = 0
    current_matching_spaces = 0
    current_deletion_1 = 0
    current_deletion_2 = 0
    current_usual_score = 0
    current_matching_letters = 0 
    current_mismatching_letters = 0
    start_gap_flag = true
    main_allignement.split("").each_with_index do |seq_1_char, position|
      seq_2_char = matching_allignement[position]
      if seq_1_char == "-" || seq_2_char == "-"
        start_gap_flag = false
        if seq_1_char == "-" && seq_2_char == "-"
          current_matching_spaces += 1
        else
          current_deletion_1 += 1 if seq_1_char == "-"
          current_deletion_2 += 1 if seq_2_char == "-"
          current_affine_score += start_gap_flag ? -4 : -2
          current_usual_score += -4
        end
      else
        start_gap_flag = true
        blosum_score = cat_cat_proxy.blosum[seq_1_char][matching_allignement[position]]
        current_affine_score += blosum_score.nil? ? 0 : blosum_score
        current_usual_score += blosum_score.nil? ? 0 : blosum_score
        if seq_1_char == matching_allignement[position]
          current_matching_letters += 1
        else
          current_mismatching_letters += 1
        end
      end
    end
    if global
      self.cat_cat_result.matching_spaces = current_matching_spaces
      self.cat_cat_result.deletions_1 = current_deletion_1
      self.cat_cat_result.deletions_2 = current_deletion_2
      self.cat_cat_result.matching_letters = current_matching_letters
      self.cat_cat_result.mismatching_letters = current_mismatching_letters
    end
    return [current_affine_score, current_usual_score]
  end

  def get_local_borders(local, seq)
    local_borders = []
    good_part = 0
    bad_part = 0
    marked = false
    seq_index = 0
    local.each_char.with_index do |char, index|
      if char == "|"
        bad_part = (bad_part+1)%2
        marked = false
        next
      end
      if char == "+"
        good_part = (good_part+1)%2
        marked = false
        local_borders << seq_index  if good_part == 0
        next
      end
      if char == "-"
        while seq[seq_index] == "-"
          seq_index += 1
        end
        next
      end
      if char != "-"
        while seq[seq_index] == "-"
          seq_index += 1
        end
        if char == seq[seq_index]
          if (!marked) && (good_part == 1)
            local_borders << seq_index
            marked = true
          end
          seq_index += 1
          next
        end
      end
    end
    # puts "L : #{local}"
    # puts "S : #{seq}"
    # puts "B : #{local_borders}"
    return local_borders
  end

  def save_to_csv(output_csv)
    self.cat_cat_result.save(output_csv, cat_cat_proxy)
  end

  def self.clear_output_file(output_filename)
    File.write("#{output_filename}_allignements.txt", '')
    File.write("#{output_filename}_local_allignements.txt", '')
  end
  
  def self.csv_header(output_csv)
    CSV.open("#{output_csv}.csv", "w") do |csv|
      header = ["pair_id",
              "Spec1",
              "Spec2",
              "exons_org_1",
              "exons_org_2",
              "Start",
              "End",
              "Leng",
              "Ex1_Sc",
              "Ex2_Sc",
              "Mono_sc",
              "Aff_sc",
              "Aff_Seq1",
              "Aff_Seq2",
              "Min_aff_seq_score",
              "Align1",
              "Align2",
              "Normed_aff_score"
            ]
      header += [ "local_score",
                  "local_aff_diff",
                  "local_score_1",
                  "local_score_2",
                  "local_score_coef",
                  "local_score_1_coef",
                  "local_score_2_coef",
                  "align_1",
                  "align_2"]
      combined_values = [
          "local_min",
          "RLoc_1",
          "RLoc_2",
          "Rmin",
          "Rmax",
      ]
      header += combined_values
      csv << header
    end
  end


end