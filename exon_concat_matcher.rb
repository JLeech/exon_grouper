require_relative "exon.rb"
require_relative "data_processor.rb"
require_relative "local_aligner.rb"

class ExonConcatMatcher

    attr_accessor :blosum
    attr_accessor :seq_1
    attr_accessor :seq_2
    attr_accessor :coords_1
    attr_accessor :coords_2

    attr_accessor :rloc_1
    attr_accessor :rloc_2

    attr_accessor :local_data

    attr_accessor :sequences
    attr_accessor :coords

    attr_accessor :added_spaces
    attr_accessor :matching_spaces
    attr_accessor :matching_letters
    attr_accessor :mismatching_letters
    
    attr_accessor :deletions_1
    attr_accessor :deletions_2

    attr_accessor :affine_score
    attr_accessor :usual_score
    attr_accessor :seq_1_score
    attr_accessor :seq_2_score

    attr_accessor :exon_1
    attr_accessor :exon_2
    attr_accessor :sequence_data

    attr_accessor :organism
    attr_accessor :match_organism
    attr_accessor :exon
    attr_accessor :match_exon

    attr_accessor :nExCat
    attr_accessor :firstEx2
    attr_accessor :lastEx2
    attr_accessor :concat_data 

    def initialize(sequences = [], coords = [], sequence_data = {}, blosum = nil, organism, match_organism, exon, match_exon)
        self.blosum = blosum.nil? ? DataProcessor.parse_blossum : blosum
        self.sequences = sequences
        self.coords = coords

        self.added_spaces = 0
        self.matching_spaces = 0
        self.matching_letters = 0
        self.mismatching_letters = 0
        self.affine_score = 0
        self.usual_score = 0

        self.deletions_1 = 0
        self.deletions_2 = 0

        self.sequence_data = sequence_data

        # NEED TO BE OPTIMAZED
        self.organism = organism
        self.match_organism = match_organism
        self.exon = exon
        self.match_exon = match_exon
        extract_data_from_sequence_data(sequence_data)
        self.concat_data = {}
    end

    def count_everything(params = {})
        parse_seq
        margin_allignements if params[:alligned].nil?
        self.affine_score, self.usual_score = count_blosum
        self.seq_1_score, _ = count_blosum(self.seq_1, self.seq_1, false)
        self.seq_2_score, _ = count_blosum(self.seq_2, self.seq_2, false)
        self.local_data = LocalAligner.new(self.seq_1, self.seq_2, self.blosum).align
        self.rloc_1 = self.local_data["local_score"].to_f/self.seq_1_score.to_f
        self.rloc_2 = self.local_data["local_score"].to_f/self.seq_2_score.to_f
        self.concat_data = get_concat_data
    end

    def parse_seq
        if sequences.empty?
            File.readlines("./data/for_exon_matcher").each_with_index do |line, index|
                self.seq_1 = line.strip if index == 0
                self.seq_2 = line.strip if index == 1
                self.coords_1 = line.strip.split(":").map(&:to_i) if index == 2
                self.coords_2 = line.strip.split(":").map(&:to_i) if index == 3
            end
        else
            self.seq_1 = self.sequences[0]
            self.seq_2 = self.sequences[1]
            self.coords_1 = self.coords[0]
            self.coords_2 = self.coords[1]
        end

        self.exon_1 = Exon.new(self.coords_1[0], self.coords_1[1], self.seq_1)
        self.exon_2 = Exon.new(self.coords_2[0], self.coords_2[1], self.seq_2)
    end

    def margin_allignements
        start_diff = coords_1[0] - coords_2[0]
        if start_diff > 0
            self.added_spaces += start_diff
            self.seq_1 = "-"*start_diff + self.seq_1
        else
            self.added_spaces += (-start_diff)
            self.seq_2 = "-"*(-start_diff) + self.seq_2
        end
        max_length = [self.seq_1.length, self.seq_2.length].max
        if self.seq_1.length < max_length
            self.added_spaces += (max_length - self.seq_1.length)
            self.seq_1 = self.seq_1 + "-"*(max_length - self.seq_1.length)
        end
        if self.seq_2.length < max_length
            self.added_spaces += (max_length - self.seq_2.length)
            self.seq_2 = self.seq_2 + "-"*(max_length - self.seq_2.length)
        end
    end

    def count_blosum(main_allignement = self.seq_1, matching_allignement = self.seq_2, global = true)
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
                blosum_score = blosum[seq_1_char][matching_allignement[position]]
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
            self.matching_spaces = current_matching_spaces
            self.deletions_1 = current_deletion_1
            self.deletions_2 = current_deletion_2
            self.matching_letters = current_matching_letters
            self.mismatching_letters = current_mismatching_letters
        end
    
        return [current_affine_score, current_usual_score]
    end

    def print_statistics_for_txt(output_filename)
        output = ""
        output += "#{self.sequence_data[:pair_id]} #{self.sequence_data[:exon_index]} #{self.sequence_data[:match_exon_index]} AFF: #{self.affine_score} Rloc1: #{self.rloc_1} Rloc2: #{self.rloc_2}\n"
        #output += "\"#{self.seq_1}\"\n"
        #output += "\"#{self.seq_2}\"\n"
        output += close_exons
        File.open("#{output_filename}_allignements.txt", 'a') { |file| file.write(output) }

        output = ""
        output += "#{self.sequence_data[:pair_id]} #{self.sequence_data[:exon_index]} #{self.sequence_data[:match_exon_index]} \n"
        output += "\"#{self.local_data['align_1']}\"\n"
        output += "\"#{self.local_data['align_2']}\"\n"
        File.open("#{output_filename}_local_allignements.txt", 'a') { |file| file.write(output) }
    end

    def close_exons
        seq_1_updated = ""
        seq_2_updated = ""

        if self.coords_1[0] > self.coords_2[0]
            seq_1_updated += organism.allignement[self.coords_2[0]..(self.coords_1[0]-1)].downcase
        elsif self.coords_1[0] < self.coords_2[0]
            seq_2_updated += match_organism.allignement[self.coords_1[0]..(self.coords_2[0]-1)].downcase
        end

        seq_1_updated += organism.allignement[self.coords_1[0]..self.coords_1[1]]
        seq_2_updated += match_organism.allignement[self.coords_2[0]..self.coords_2[1]]

        if self.coords_1[1] > self.coords_2[1]
            seq_2_updated += match_organism.allignement[(self.coords_2[1]+1)..self.coords_1[1]].downcase
        elsif self.coords_1[1] < self.coords_2[1]
            seq_1_updated += organism.allignement[(self.coords_1[1]+1)..self.coords_2[1]].downcase
        end                

        return "\"#{seq_1_updated}\"\n\"#{seq_2_updated}\"\n"

        #puts organism.allignement[self.coords_1[0]..self.coords_1[1]]
        #puts self.seq_2
    end

    def extract_data_from_sequence_data(sequence_data)
        concatted_exons = sequence_data[:concatted_exons]
        self.nExCat = concatted_exons.length
        self.firstEx2 = concatted_exons.first.uuid
        self.lastEx2 = concatted_exons.last.uuid
    end

  def get_concat_data
      counter = 0
      data = {"counter"=> 0,
              "first_domain_index"=> "", 
              "last_domain_index" => "",
              "first_coverage_per" => "",
              "last_coverage_per" => ""
            }
      parts = self.local_data["align_2"].split("UU").map { |part| part.split("uu") }.flatten
      vals = parts.map { |part| part.chars.map(&:ord).map { |place| place <= 90  } } #90 = "Z"
      has_uppercase = vals.map { |part| part.uniq.include?(true) }
      data["counter"] = has_uppercase.count(true)
      return data if data["counter"] == 0
      first_domain_pre_index = has_uppercase.index(true)
      data["first_domain_index"] = self.sequence_data[:concatted_exons][first_domain_pre_index].uuid
      last_domain_pre_index = has_uppercase.length - has_uppercase.reverse.index(true)-1
      data["last_domain_index"] = self.sequence_data[:concatted_exons][last_domain_pre_index].uuid

      upper_count_first = vals[first_domain_pre_index].count(true)
      data["first_coverage_per"] = upper_count_first/vals[first_domain_pre_index].length.to_f*100
      upper_count_last = vals[last_domain_pre_index].count(true)
      data["last_coverage_per"] = upper_count_last/vals[last_domain_pre_index].length.to_f*100
      return data
  end

    def print_for_csv(output_filename)
        exon_1_length = self.exon_1.allignement.gsub("-","").length
        exon_2_length = self.exon_2.allignement.gsub("-","").length
        data = [
                self.sequence_data[:pair_id],
                self.sequence_data[:org_name],
                self.sequence_data[:exon_index],
                self.coords_1[0],
                self.coords_1[1],
                exon_1_length,
                self.seq_1_score.to_f.round(2),
                self.sequence_data[:match_org_name],
                self.sequence_data[:match_exon_index],
                self.coords_2[0],
                self.coords_2[1],
                exon_2_length,
                self.nExCat,
                self.firstEx2,
                self.lastEx2,
                self.concat_data["counter"],
                self.concat_data["first_domain_index"],
                self.concat_data["last_domain_index"],
                self.concat_data["first_coverage_per"],
                self.concat_data["last_coverage_per"],
                self.seq_2_score.to_f.round(2),
                self.affine_score,
                self.usual_score,
                self.added_spaces,
                self.matching_spaces,
                self.matching_letters,
                self.mismatching_letters,
                self.deletions_1,
                self.deletions_2,
                self.exon_1.get_exons_matching_coords(self.exon_2)-1,
                (self.affine_score/self.seq_1_score.to_f).round(2),
                (self.affine_score/self.seq_2_score.to_f).round(2),
                [self.affine_score/self.seq_1_score.to_f,self.affine_score/self.seq_2_score.to_f].min.round(2),
                ([exon_1_length,exon_2_length].min/[exon_1_length,exon_2_length].max.to_f).round(2),
                (self.matching_letters / exon_1_length.to_f).round(2),
                (self.matching_letters / exon_2_length.to_f).round(2),
                "",
                "",
                (self.affine_score/[self.seq_2_score,self.seq_1_score].max.to_f).round(2)
            ]
        data += self.local_data.values
        data += [
            self.rloc_1,
            self.rloc_2,
            [self.rloc_1,self.rloc_2].min,
            [self.rloc_1,self.rloc_2].max
        ]
        CSV.open("#{output_filename}.csv", "a") { |csv| csv << data}
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
                  "ConcatEx",
                  "Start2",
                  "End2",
                  "Leng2",
                  "NExCat",
                  "1stEx2",
                  "LastEx2",
                  "Dom2NEx",
                  "1stExDom",
                  "LastExDom",
                  "Cover1st",
                  "CoverLast",
                  "ConcatEx_Sc",
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
                        "align_2",
                        "local_min",]
            combined_values = [
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