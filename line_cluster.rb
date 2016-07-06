class LineCluster

	attr_accessor :starts
	attr_accessor :ends
	attr_accessor :length
	attr_accessor :clusters

	def initialize(borders)
		self.starts = borders.map(&:first)
		self.ends = borders.map(&:last)
		self.lengths = [ends,starts].transpose.map {|x| x.reduce(:-)}
	end 

	def cluster
		
	end

end