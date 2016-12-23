# # NExCat  Количество экзонов в катенации      
# # 1stEx2  Номер 1-го экзона       
# # LastEx2 Номер последнего экзона     
# # Dom2 #Ex    Количество экзонов в проекции BLS на 2-й организм       
# # 1stExDom    Номер 1-го экзона в этой проекции       
# # LfstExDom   Номер последнего экзона     
# # Cover1st    Доля покрытия 1-го экзона       
# # CoverLast   Доля покрытия последнего экзона     

# def get_concat_data(str)
#     counter = 0
#     data = {"counter"=> 0,
#             "first_domain_index"=> "", 
#             "last_domain_index" => "",
#             "first_coverage_per" => "",
#             "last_coverage_per" => ""
#           }
#     parts = str.split("UU").map { |part| part.split("uu") }.flatten
#     vals = parts.map { |part| part.chars.map(&:ord).map { |place| place <= 90  } }
#     has_uppercase = vals.map { |part| part.uniq.include?(true) }
#     data["counter"] = has_uppercase.count(true)
#     return "#{data}" if data["counter"] == 0
#     data["first_domain_index"] = has_uppercase.index(true)

#     data["last_domain_index"] = has_uppercase.length - has_uppercase.reverse.index(true)-1
#     upper_count_first = vals[data["first_domain_index"]].count(true)
#     data["first_coverage_per"] = upper_count_first/vals[data["first_domain_index"]].length.to_f*100
#     upper_count_last = vals[data["last_domain_index"]].count(true)
#     data["last_coverage_per"] = upper_count_first/vals[data["last_domain_index"]].length.to_f*100
#     return "#{data}"
# end
# x = "xxZ"
# puts x
# puts get_concat_data(x)
