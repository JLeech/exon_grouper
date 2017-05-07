require_relative "exon.rb"
require_relative "data_processor.rb"
require_relative "cat_cat_local_aligner.rb"
require_relative "common.rb"


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
    local_1_for_score = ""
    local_2_for_score = ""
    locals.each do |local_split|
      local_1 += local_split.get_for_coords_seq(1)
      local_2 += local_split.get_for_coords_seq(2)
    end
    locals.each do |local_split|
      next if local_split.type == LocalReqursiveResult::BAD 
      local_1_for_score += local_split.get_raw.first
      local_2_for_score += local_split.get_raw.last
    end

    locals.each do |local_split|
      next if local_split.type == LocalReqursiveResult::BAD
      local_split.score, _ = count_blosum(local_split.seq_1,local_split.seq_2)
      local_split.score_seq_1, _ = count_blosum(local_split.seq_1,local_split.seq_1)
      local_split.score_seq_2, _ = count_blosum(local_split.seq_2,local_split.seq_2)
    end

    set_exons_for_locals(locals)
    IterSaver.save(cat_cat_proxy, locals)

    cat_cat_result.locals = locals
    cat_cat_result.local_score, _ = count_blosum(local_1_for_score, local_2_for_score)
    cat_cat_result.local_self_1_score, _ = count_blosum(local_1_for_score, local_1_for_score)
    cat_cat_result.local_self_2_score, _ = count_blosum(local_2_for_score, local_2_for_score)
    cat_cat_result.local_borders_seq_1 = get_local_borders(local_1, cat_cat_proxy.seq_1)
    cat_cat_result.local_borders_seq_2 = get_local_borders(local_2, cat_cat_proxy.seq_2)
    raise "local 1 borders error" if cat_cat_result.local_borders_seq_1.length%2 != 0
    raise "local 2 borders error" if cat_cat_result.local_borders_seq_2.length%2 != 0
  end

  def has_multiple?
    return (cat_cat_result.locals.map(&:seq1_exons_ids).map(&:length).max > 1 ) & (cat_cat_result.locals.map(&:seq2_exons_ids).map(&:length).max > 1)
  end

  def local_recursive(seq_1, seq_2, len_coef = 1.0, iter = 1)
    if stop_results?(seq_1, seq_2, len_coef)
      al_1, al_2 = margin(seq_1, seq_2)
      return [LocalReqursiveResult.new(al_1, al_2, LocalReqursiveResult::BAD)]
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
      al_1, al_2 = margin(seq_1, seq_2)
      final = [LocalReqursiveResult.new(al_1, al_2, LocalReqursiveResult::BAD)]
    else
      self.cat_cat_result.local_iters += 1
      final = local_recursive(left_part_1,left_part_2,left_length_coef, iter+1) + 
              [LocalReqursiveResult.new(local_result.align_1, local_result.align_2, LocalReqursiveResult::GOOD, iter)] + 
              local_recursive(right_part_1,right_part_2,right_length_coef,iter+1)
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
    return [seq_1, seq_2]
  end

  def set_exons_for_locals(locals)
    cur_org_ex_index = 0 
    cur_local_index = 0

    seq_1 = locals.map(&:seq_1).join("_")

    marked = false
    seq_1.each_char.with_index do |char, index|
      if char == "_"
        cur_local_index += 1
        marked = false
        next
      end
      if (seq_1[index..(index+1)] == "UU") & (index != 0)  
        cur_org_ex_index += 1
        marked = false
        next
      end
      if char != "U" && marked == false
        locals[cur_local_index].seq1_exons_ids << cat_cat_proxy.get_org_exons_raw[cur_org_ex_index]
        marked = true
      end
    end

    cur_org_ex_index = 0 
    cur_local_index = 0

    seq_2 = locals.map(&:seq_2).join("_")

    marked = false
    seq_2.each_char.with_index do |char, index|
      if char == "_"
        cur_local_index += 1
        marked = false
        next
      end
      if (seq_2[index..(index+1)] == "UU") & (index != 0)  
        cur_org_ex_index += 1
        marked = false
        next
      end
      if char != "U" && marked == false
        locals[cur_local_index].seq2_exons_ids << cat_cat_proxy.get_match_org_exons_raw[cur_org_ex_index]
        marked = true
      end
    end

  end

  def count_blosum(main_alignment, matching_alignment, global = false)
    current_affine_score = 0
    current_matching_spaces = 0
    current_deletion_1 = 0
    current_deletion_2 = 0
    current_usual_score = 0
    current_matching_letters = 0 
    current_mismatching_letters = 0
    start_gap_flag = true
    main_alignment.split("").each_with_index do |seq_1_char, position|
      seq_2_char = matching_alignment[position]
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
        blosum_score = cat_cat_proxy.blosum[seq_1_char][matching_alignment[position]]
        current_affine_score += blosum_score.nil? ? 0 : blosum_score
        current_usual_score += blosum_score.nil? ? 0 : blosum_score
        if seq_1_char == matching_alignment[position]
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
        local_borders << seq_index-1  if good_part == 0
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
    return local_borders
  end

  def save_to_csv(output_csv)
    self.cat_cat_result.save(output_csv, cat_cat_proxy)
  end

  def self.clear_output_file(output_filename)
    File.write("#{output_filename}_alignments.txt", '')
    File.write("#{output_filename}_local_alignments.txt", '')
  end


  def self.csv_header(output_csv)
    CSV.open("#{output_csv}.csv", "w") do |csv|
      header = ["pair_id",
              "Spec1",
              "Spec2",
              "Spec1_code_name",
              "Spec2_code_name",
              "exons_org_1",
              "exons_org_2",
              "exons_org_1_count",
              "exons_org_2_count",
              "maxS",
              "minS",
              "Start",
              "End",
              "Leng",
              "Tuple1_Sc",
              "Tuple2_Sc",
              "Aff_sc",
              "Aff_Seq1",
              "Aff_Seq2",
              "Min_aff_seq_score",
              "Align1",
              "Align2"
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