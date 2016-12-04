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

  def initialize(sequences = [], coords = [], blosum = nil, organism, match_organism, pair_id)
    self.sequences = sequences
    self.coords = coords
    self.blosum = blosum.nil? ? DataProcessor.parse_blossum : blosum
    self.organism = organism
    self.match_organism = match_organism
    self.pair_id = pair_id
  end

  def seq_1
    self.sequences[0]
  end

  def seq_2 
    self.sequences[1]
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

  attr_accessor :local_data

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
  end

  def check
    if ((seq_1_score < usual_score) || (seq_2_score < usual_score))
      raise "self score > usual score"  
    end

  end

  def print
    puts "#{added_spaces}"
    puts "#{matching_spaces}"
    puts "#{matching_letters}"
    puts "#{mismatching_letters}"
    puts "#{affine_score}"
    puts "#{usual_score}"
    puts "#{deletions_1}"
    puts "#{deletions_2}"
    puts "#{seq_1_score}"
    puts "#{seq_2_score}"
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
    #cat_cat_result.local_data = LocalAligner.new(cat_cat_proxy.seq_1, cat_cat_proxy.seq_2, cat_cat_proxy.blosum).align
    make_local_recursive
    cat_cat_result.check
  end

  def make_local_recursive
    no_gap_seq_1 = cat_cat_proxy.seq_1.gsub("-","")
    no_gap_seq_2 = cat_cat_proxy.seq_2.gsub("-","")
    local_data = CatCatLocalAligner.new(no_gap_seq_1, no_gap_seq_2, cat_cat_proxy.blosum).align
    puts "#{cat_cat_proxy.seq_1}"
    puts "#{cat_cat_proxy.seq_2}"
    puts "#{local_data['align_1']}"
    puts "#{local_data['align_2']}"
    puts "--------------------------"
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
          current_affine_score += start_gap_flag ? -9 : -2
          current_usual_score += -9
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

  def self.clear_output_file(output_filename)
    File.write("#{output_filename}_allignements.txt", '')
    File.write("#{output_filename}_local_allignements.txt", '')
  end
  
  def self.csv_header(output_csv)
    CSV.open("#{output_csv}.csv", "w") do |csv|
      header = ["pair_id",
              "Spec1",
              "Ex1",
              "Start1",
              "End1",
              "Leng1",
              "Ex1_Sc",
              "Spec2",
              "Ex2",
              "Start2",
              "End2",
              "Leng2",
              "Ex2_Sc",
              "Aff_sc",
              "Mono_sc",
              "Add_sp",
              "Match_sp",
              "Match",
              "MisM",
              "Del1",
              "Del2",
              "Inters",
              "Aff_Seq1",
              "Aff_Seq2",
              "Min_aff_seq_score",
              "Leng_Min_Max",
              "Per_id1",
              "Per_id2",
              "Align1",
              "Align2",
              "Normed_aff_score"
            ]
      header += [ "start_position_1",
                  "end_position_1",
                  "right_margin_1",
                  "start_position_2",
                  "end_position_2",
                  "right_margin_2",
                  "local_score",
                  "local_score_1",
                  "local_score_2",
                  "local_score_coef",
                  "local_length_coef",
                  "local_length_1_coef",
                  "local_length_2_coef",
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