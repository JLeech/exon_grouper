class ConCatMatcher

  attr_accessor :proxy

  def initialize(proxy)
    self.proxy = proxy
  end


  def align
    
  end

  def save
  
  end

  def self.header(outpu_filname)
    header = [
              "pair_id",
              "Spec1_code_name",
              "Spec2_code_name",
              "exons_org_1",
              "exons_org_2",
              "current_exon",
              "current_exons_match",
              "first_domain_index",
              "last_domain_index",
              "first_cover_percent",
              "last_cover_percent",
              "local_score",
              "self_local_score_1",
              "self_local_score_2",
              "local_score_1_coef",
              "local_score_2_coef",
    ]
  end




end