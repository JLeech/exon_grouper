require_relative "exon.rb"
require_relative "data_processor.rb"
require_relative "local_aligner.rb"

class ExonMatcher

    attr_accessor :blosum
    attr_accessor :seq_1
    attr_accessor :seq_2
    attr_accessor :coords_1
    attr_accessor :coords_2

    # attr_accessor :local_score
    # attr_accessor :local_seq_1_score
    # attr_accessor :local_seq_2_score
    # attr_accessor :local_pair

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

    def initialize(sequences = [], coords = [], sequence_data = {}, blosum = nil)
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
    end

    def count_everything(params = {})
        parse_seq
        margin_allignements if params[:alligned].nil?
        self.affine_score, self.usual_score = count_blosum
        self.seq_1_score, _ = count_blosum(self.seq_1, self.seq_1, false)
        self.seq_2_score, _ = count_blosum(self.seq_2, self.seq_2, false)
        self.local_data = LocalAligner.new(self.seq_1, self.seq_2, self.blosum).align
        #self.local_score, self.local_seq_1_score, self.local_seq_2_score, self.local_pair  = ExonMatcher.get_local_score(self.seq_1, self.seq_2)
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
        output += "#{self.sequence_data[:pair_id]}\n"
        output += "\"#{self.seq_1}\"\n"
        output += "\"#{self.seq_2}\"\n"
        File.open("#{output_filename}_allignements.txt", 'a') { |file| file.write(output) }

        output = ""
        output += "#{self.sequence_data[:pair_id]}\n"
        output += "\"#{self.local_data['align_1']}\"\n"
        output += "\"#{self.local_data['align_2']}\"\n"
        File.open("#{output_filename}_local_allignements.txt", 'a') { |file| file.write(output) }
    end

    def print_for_csv(output_filename)
        exon_1_length = self.exon_1.allignement.length
        exon_2_length = self.exon_2.allignement.length

        #local_leng_1 = get_no_gap_langth(local_seq_1)
        #local_leng_2 = get_no_gap_langth(local_seq_2)

        #local_aff_score = count_blosum(local_seq_1,local_seq_2,false)[0]
        #loc_seq_1_score = count_blosum(local_seq_1,local_seq_1,false)[0]
        #loc_seq_2_score = count_blosum(local_seq_2,local_seq_2,false)[0]
        #aff_loc_seq_1_score = (local_aff_score/loc_seq_1_score.to_f).round(2)
        #aff_loc_seq_2_score = (local_aff_score/loc_seq_2_score.to_f).round(2)
        #if "#{aff_loc_seq_1_score}" == "NaN"
        #  loc_seq_1_score = -999999
        #  loc_seq_2_score = -999999
        #  aff_loc_seq_1_score = -999999
        #  aff_loc_seq_2_score = -999999
        #end
        data = [self.sequence_data[:pair_id],
                self.sequence_data[:org_name],
                self.sequence_data[:exon_index] + 1,
                self.coords_1[0],
                self.coords_1[1],
                self.coords_1.reverse.inject(:-),
                self.seq_1_score.to_f.round(2),
                self.sequence_data[:match_org_name],
                self.sequence_data[:match_exon_index] + 1,
                self.coords_2[0],
                self.coords_2[1],
                self.coords_2.reverse.inject(:-),
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
                # self.local_score,
                # self.local_seq_1_score,
                # self.local_seq_2_score,
                # self.local_pair,
                (self.affine_score/[self.seq_2_score,self.seq_1_score].max.to_f).round(2)
                # local_leng_1,
                # local_leng_2,
                # ([local_leng_1,local_leng_2].min/[local_leng_1,local_leng_2].max.to_f).round(2),                
                # local_aff_score,
                # loc_seq_1_score,
                # loc_seq_2_score,
                # aff_loc_seq_1_score,
                # aff_loc_seq_2_score,
                # [aff_loc_seq_1_score, aff_loc_seq_2_score].min
            ]
        data += self.local_data.values
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
                  # "Local_score",
                  # "local_seq_1_score",
                  # "local_seq_2_score",
                  # "local_pair",
                  "Normed_aff_score"
                  # "Leng_Al",
                  # "Local_leng_1",
                  # "Local_leng_2",
                  # "Local_leng_ratio",
                  # "Local_aff_score",
                  # "Loc_seq_1_score",
                  # "Loc_seq_2_score",
                  # "Aff_loc_seq_1_score",
                  # "Aff_loc_seq_2_score",
                  # "min_Aff_loc_score"
                ]
            header += [ "start_position_1",
                        "start_position_2",
                        "end_position_1",
                        "end_position_2",
                        "score_1",
                        "score_2",
                        "local_length_coef",
                        "local_score_1_coef",
                        "local_score_2_coef"]
            csv << header
        end
    end

##### local block
    # def self.get_local_score(sequence_1, sequence_2)
    #   score_all = count_local_score_for_sequences(sequence_1, sequence_2)
    #   score_1 = count_local_score_for_sequences(sequence_1, sequence_1)
    #   score_2 = count_local_score_for_sequences(sequence_2, sequence_2)
    #   return [ score_all/[score_1, score_2].max, score_1, score_2, score_all] 
    # end

    # def self.count_local_score_for_sequences(sequence_1, sequence_2)
    #   path_1 = "data/1.fasta"
    #   path_2 = "data/2.fasta"
    #   File.write(path_1,">1\n#{sequence_1.gsub('-','')}")
    #   File.write(path_2,">2\n#{sequence_2.gsub('-','')}")
    #   return `./ssw_test -p #{path_1} #{path_2}`.to_f
    # end

    # def get_local_allignement(sequence_1, sequence_2)
    #     chars_1,chars_2 = form_islands(sequence_1, sequence_2)

    #     chars_1,chars_2 = delete_till_both_letters( chars_1.reverse, chars_2.reverse )
    #     #return [ chars_1.reverse.join, chars_2.reverse.join ]
    # end



#####

private

    # def form_islands(sequence_1, sequence_2)
    #   chars_1, chars_2 = delete_till_both_letters( sequence_1.split(""), sequence_2.split("") )
    #   place = get_island_end(chars_1, chars_2)
    #   unless place == chars_1.length
    #     island_score = count_blosum(chars_1[0..place].join(""),chars_2[0..place].join(""))[0]
    #     river_end = get_river_end(chars_1, chars_2, place)
    #     river_score = count_blosum(chars_1[place..river_end].join(""),chars_2[place..river_end].join(""))[0]
    #     if 
    #   end
    # end

    # def get_no_gap_langth(string)
    #   return string.gsub("-","").length
    # end

    # def get_river_end(chars_1, chars_2, place)
    #   new_chars_1 = chars_1[0..place]
    #   new_chars_2 = chars_2[0..place]
    #   new_chars_2.length.times do |river_place|
    #     if new_chars_1[river_place] != "-" && new_chars_2[river_place] != "-"
    #       return place + river_place
    #     end
    #   end
    #   return chars_1.length
    # end

    # def get_island_end(chars_1, chars_2)
    #   chars_1.length.times do |place|
    #     if chars_1[place] == "-" && chars_2[place] == "-"
    #       return place
    #     end
    #   end
    #   return chars_1.length
    # end

    # def delete_till_both_letters(chars_1, chars_2)
    #     while true
    #         if chars_1.last == "-" || chars_2.last == "-"
    #             chars_1.pop
    #             chars_2.pop
    #         else
    #             break
    #         end
    #         break if chars_1.empty? || chars_2.empty?
    #     end
    #     return [chars_1, chars_2]
    # end

end
