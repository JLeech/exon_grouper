require_relative "exon.rb"

class Organism

  attr_accessor :name
  attr_accessor :code_name
  attr_accessor :exons
  attr_accessor :exons_coords
  attr_accessor :number
  attr_accessor :alignment_length
  attr_accessor :alignment
  attr_accessor :points

  def initialize(name, exons, number, alignment)
    self.name = name
    self.exons = exons
    self.number = number
    self.alignment_length = alignment.length
    self.alignment = alignment
    self.code_name = name.split("_").map { |val| val[0..2] }.join("_")
    self.exons_coords = set_exons_hash
  end

  def save_references(file_path)
    data = ""
    data += "#{self.name},#{self.code_name}, #{self.exons.length}\n"
    File.open("#{file_path}_references.csv", 'a') { |file| file.write(data) }
  end

  def set_exons_hash
    # {exon_uid -> [exon_coords]}
    ex_hash =  {}
    exons.each { |exon| ex_hash[exon.uuid.to_s] = exon.get_coords }
    return ex_hash
  end

  def self.set_headers(file_path)
    header = "organism name, organism code, #exons\n"
    File.open("#{file_path}_references.csv", 'w') { |file| file.write(header) }
  end
end





