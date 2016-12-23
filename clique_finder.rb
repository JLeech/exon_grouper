# require 'set'

# class Set
#   def take!(args)
#     taken = self.take(args)
#     self.subtract(taken)
#     return taken
#   end
# end

# class Kerbosh

#     attr_accessor :cliques
#     attr_accessor :graph

#     def initialize(graph)
#       self.cliques = []
#       self.graph = graph
#     end

#     def find_cliques
#       p = Set.new(graph.keys)
#       r = Set.new()
#       x = Set.new()
#       degeneracy_ordering(graph).each do |vert|
#         neighs = graph[vert]
#         find_cliques_pivot(graph, r.union([vert]), p.intersection(neighs), x.intersection(neighs))
#         p -= [vert]
#         x += [vert]
#       end
#       return self.cliques.map(&:to_a)
#     end

#     def degeneracy_ordering(graph)
#       ordering = []
#       ordering_set = Set.new()
#       degrees = Hash.new(0)
#       degen = Hash.new {|hsh, key| hsh[key] = [] }
#       max_deg = -1
#       graph.keys.each do |vert|
#         deg = graph[vert].length
#         degen[deg] << vert
#         degrees[vert] = deg
#         if deg > max_deg
#           max_deg = deg
#         end
#       end
#       while true
#         i = 0
#         while i <= max_deg
#           if degen[i].length != 0
#             break
#           end
#           i += 1
#         end
#         break if i > max_deg
#         v = degen[i].pop
#         ordering << v
#         ordering_set.add(v)
#         graph[v].each do |w|
#           if !ordering_set.include?(w)
#             deg = degrees[w]
#             degen[deg].delete(w)
#             if deg > 0
#               degrees[w] -= 1
#               degen[deg-1] << w
#             end
#           end
#         end
#       end
#       return ordering.reverse
#     end

#     def find_cliques_pivot(graph, r, p, x)
#       if ((p.length == 0) & (x.length == 0))
#         self.cliques << r
#       else
#         u = p.union(x).first
#         p.difference(graph[u]).each do |v|
#           neighs = graph[v]
#           find_cliques_pivot(graph, r.union([v]), p.intersection(neighs), x.intersection(neighs))
#           p.delete(v)
#           x.add(v)
#         end
#       end
#     end
# end


# class CliqueFinder

#   attr_accessor :graph
#   attr_accessor :cliques
#   attr_accessor :prime_exon

#   def initialize(exon)
#     self.graph = {}
#     self.prime_exon = exon
#   end

#   def make_graph
#     add_exon_to_graph(self.prime_exon)  
#   end

#   def make_cliques
#     ker = Kerbosh.new(self.graph)
#     self.cliques = ker.find_cliques
#   end

# private

#   def add_exon_to_graph(exon)
#     graph[exon.uuid] = exon.connections.select(&:green?).map(&:uuid)
#     exon.connections.select(&:green?).each do |connect|
#       add_exon_to_graph(connect) unless graph.keys.include?(connect.uuid)
#     end
#   end

#   def print_graph
#     puts "#{graph}"
#   end

# end