bhattgenes= DataFrames.DataFrame(
    "GeneID" => Bhatt2012.inducible_genes_figure3()[!,:GeneSymbol], 
    "Class"=>Bhatt2012.inducible_genes_figure3()[!,:Class]
    );


# Bool for 2H LPS
bool1 = sce.rowData[!,:WT_2H__alpha].>=alpha
# Bool for 8H LPS
bool2 = sce.rowData[!,:WT_8H__alpha].>=alpha
# Bool for either in 2H or 8H
bool = (bool1.+bool2).>0


subsubsce = filter_genes(bool,sce)

sceBhatt = innerjoin(subsubsce.rowData, bhattgenes,on = :GeneID)
sceBhatt[!,:Class] = [replace(replace(ii, "A1"=>"A1+2"), "A2"=>"A1+2") for ii in  sceBhatt[!,:Class]]

println(string("Percent ", alpha*100, "%"))
println(string("Total genes considered ", nrow(sceBhatt)))


p1s = []
p1as = []

for class in sort!(unique(sceBhatt[!,"Class"]))
    #println("Class $class")
    #println("ln CPM plus 1 in expressing cells")
   
    #println("WT UT vs WT 2H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_UT__mu"]], [ii for ii in sub[!,"WT_2H__mu"]])
    #println(t)
    p1 = pvalue(t)  
    push!(p1s, p1)
    
      #println("WT 2H vs WT 8H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_2H__mu"]], [ii for ii in sub[!,"WT_8H__mu"]])
    #println(t)
    p1 = pvalue(t)  
    push!(p1s, p1)
    
    #println("WT UT vs WT 8H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_UT__mu"]], [ii for ii in sub[!,"WT_8H__mu"]])
    #println(t)
    p1 = pvalue(t)  
    push!(p1s, p1)
    
  
    
    #println("Fraction Expressing cells")
  # println("WT UT vs WT 2H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_UT__alpha"]], [ii for ii in sub[!,"WT_2H__alpha"]])
    #println(t)
    p1a = pvalue(t)  
   push!(p1as, p1a)
    
      #println("WT 2H vs WT 8H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_2H__alpha"]], [ii for ii in sub[!,"WT_8H__alpha"]])
    #println(t)
    p1a = pvalue(t)  
    push!(p1as, p1a)
    
    #println("WT UT vs WT 8H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_UT__alpha"]], [ii for ii in sub[!,"WT_8H__alpha"]])
    #println(t)
    p1a = pvalue(t)  
    push!(p1as, p1a)
     
end

