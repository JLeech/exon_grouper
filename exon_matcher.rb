require_relative "exon.rb"
require_relative "data_processor.rb"

class ExonMatcher

    attr_accessor :blosum
    attr_accessor :seq_1
    attr_accessor :seq_2
    attr_accessor :coords_1
    attr_accessor :coords_2

    attr_accessor :sequences
    attr_accessor :coords

    attr_accessor :added_spaces
    attr_accessor :matching_spaces
    attr_accessor :matching_letters
    attr_accessor :not_mathing_letters
    
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
        self.not_mathing_letters = 0
        self.affine_score = 0
        self.usual_score = 0
        self.sequence_data = sequence_data
    end

    def count_everything
        parse_seq
        margin_allignements
        self.affine_score, self.usual_score = count_blosum

        self.seq_1_score, _ = count_blosum(self.seq_1, self.seq_1, false)

        self.seq_2_score, _ = count_blosum(self.seq_2, self.seq_2, false)

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
            self.seq_1 = "-"*start_diff + self.seq_1
            self.added_spaces += start_diff
        else
            self.seq_2 = "-"*(-start_diff) + self.seq_2
            self.added_spaces += (-start_diff)
        end
        max_length = [self.seq_1.length, self.seq_2.length].max
        if self.seq_1.length < max_length
            self.seq_1 = self.seq_1 + "-"*(max_length - self.seq_1.length)
            self.added_spaces += (max_length - self.seq_1.length)
        end
        if self.seq_2.length < max_length
            self.seq_2 = self.seq_2 + "-"*(max_length - self.seq_2.length)
            self.added_spaces += (max_length - self.seq_1.length)
        end
    end

    def count_blosum(main_allignement = self.seq_1, matching_allignement = self.seq_2, global = true)
        affine_score_line = ""
        seq_1_line = ""
        seq_2_line = ""
        start_gap_flag = true
        main_allignement.split("").each_with_index do |seq_1_char, position|
            seq_2_char = matching_allignement[position]
            if seq_1_char == "-" || seq_2_char == "-"
                if seq_1_char == "-" && seq_2_char == "-"
                    self.matching_spaces += 1 if global
                    start_gap_flag = false  
                else
                    self.affine_score += start_gap_flag ? -9 : -2
                    start_gap_flag = false
                    self.usual_score += -9
                end
            else
                start_gap_flag = true
                blosum_score = blosum[seq_1_char][matching_allignement[position]]
                self.affine_score += blosum_score.nil? ? 0 : blosum_score
                self.usual_score += blosum_score.nil? ? 0 : blosum_score
                if global
                    if seq_1_char == matching_allignement[position]
                        self.matching_letters += 1
                    else
                        self.not_mathing_letters += 1
                    end
                end
            end
            #affine_score_line += "#{affine_score.to_s}|"
            #seq_1_line += " "*(affine_score.to_s.length-1) + seq_1_char + "|"
            #seq_2_line += " "*(affine_score.to_s.length-1) + seq_2_char + "|"
        end

        
        #puts "_______________________"
        # puts affine_score_line
        # puts seq_1_line
        # puts seq_2_line
        return [affine_score, usual_score]
    end

    def print_statistics_for_txt(output_filename, aditional_data = {})
        exon_1_length = self.exon_1.allignement.length
        exon_2_length = self.exon_2.allignement.length
        output = ""
        output += "_________________________\n"
        output += "pair_id: #{sequence_data[:pair_id]}\n"

        output += "#{sequence_data[:org_name]} : exon_number:[#{sequence_data[:exon_index] + 1}]\n"
        output += "#{sequence_data[:match_org_name]} : exon_number:[#{sequence_data[:match_exon_index] + 1}]\n"
        output += "#{self.coords_1}\n"
        output += "#{self.coords_2}\n"

        output += self.seq_1 + "\n"
        output += self.seq_2 + "\n"

        output += "affine_score:         #{self.affine_score}\n"
        output += "usual_score:      #{self.usual_score}\n"
        output += "added_spaces:         #{self.added_spaces}\n"
        output += "matching_spaces:  #{self.matching_spaces}\n"
        output += "matching_letters:     #{self.matching_letters}\n"
        output += "not_mathing_letters:  #{self.not_mathing_letters}\n"
        output += "matching_coords:      #{self.exon_1.get_exons_matching_coords(self.exon_2)}\n"

        output += "aff/seq1_score:   #{self.affine_score/seq_1_score.to_f}\n"
        output += "aff/seq2_score:   #{self.affine_score/self.seq_2_score.to_f}\n"     

        output += "length coef:  #{[exon_1_length,exon_2_length].min.to_f/[exon_1_length,exon_2_length].max}\n"
        output += "match_letter_1 :  #{ self.matching_letters.to_f / exon_1_length } \n"
        output += "match_letter_2 :  #{ self.matching_letters.to_f / exon_2_length } \n"
        

        output += "_________________________\n\n"
        File.open("#{output_filename}.txt", 'a') { |file| file.write(output) }
    end

    def print_for_csv(output_filename)
        exon_1_length = self.exon_1.allignement.length
        exon_2_length = self.exon_2.allignement.length
        data = [self.sequence_data[:pair_id],
                self.affine_score,
                self.usual_score,
                self.added_spaces,
                self.matching_spaces,
                self.matching_letters,
                self.not_mathing_letters,
                self.exon_1.get_exons_matching_coords(self.exon_2),
                self.affine_score/seq_1_score.to_f,
                self.affine_score/self.seq_2_score.to_f,
                [exon_1_length,exon_2_length].min.to_f/[exon_1_length,exon_2_length].max,
                self.matching_letters / exon_1_length.to_f,
                self.matching_letters / exon_2_length.to_f]
        CSV.open("#{output_filename}.csv", "a") { |csv| csv << data}
    end
    
    def self.clear_output_file(output_filename)
        File.write("#{output_filename}.txt", '')
    end

    def self.csv_header(output_csv)
        CSV.open("#{output_csv}.csv", "w") do |csv|
          csv << ["pair_id",
                  "affine_score", 
                  "usual_score", 
                  "added_spaces", 
                  "matching_spaces", 
                  "matching_letters", 
                  "not_mathing_letters", 
                  "matching_coords", 
                  "aff/seq1_score", 
                  "aff/seq2_score", 
                  "length_coef", 
                  "match_letter_1", 
                  "match_letter_2"]
        end
    end
        

    def print_allignements
        puts "#{seq_1}"
        puts "#{seq_2}"
        puts "_______________________"
    end

    def print_coords
        puts "seq_1 coords: #{coords_1}"
        puts "seq_2 coords: #{coords_2}"
        puts "_______________________"
    end

end




# exon_matcher = ExonMatcher.new
# exon_matcher.count_everything

#exon_matcher.print_coords
#exon_matcher.print_allignements
#exon_matcher.margin_allignements
#exon_matcher.print_allignements

#exon_matcher.count_blosum_and_print
#exon_matcher.print_statistics

