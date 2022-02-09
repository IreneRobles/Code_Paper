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
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_UT__bm"]], [ii for ii in sub[!,"WT_2H__bm"]])
    #println(t)
    p1 = pvalue(t)  
    push!(p1s, p1)
    
      #println("WT 2H vs WT 8H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_2H__bm"]], [ii for ii in sub[!,"WT_8H__bm"]])
    #println(t)
    p1 = pvalue(t)  
    push!(p1s, p1)
    
    #println("WT UT vs WT 8H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_UT__bm"]], [ii for ii in sub[!,"WT_8H__bm"]])
    #println(t)
    p1 = pvalue(t)  
    push!(p1s, p1)
    
  
    
    #println("Fraction Expressing cells")
  # println("WT UT vs WT 2H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_UT__fm"]], [ii for ii in sub[!,"WT_2H__fm"]])
    #println(t)
    p1a = pvalue(t)  
   push!(p1as, p1a)
    
      #println("WT 2H vs WT 8H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_2H__fm"]], [ii for ii in sub[!,"WT_8H__fm"]])
    #println(t)
    p1a = pvalue(t)  
    push!(p1as, p1a)
    
    #println("WT UT vs WT 8H")
    sub = sceBhatt[sceBhatt[!,"Class"].==class, :]
    t = HypothesisTests.SignedRankTest([ii for ii in sub[!,"WT_UT__fm"]], [ii for ii in sub[!,"WT_8H__fm"]])
    #println(t)
    p1a = pvalue(t)  
    push!(p1as, p1a)
     
end


figure(figsize = (20, 10))


subplot(1,2,1)
y = "moment burst frequency"
pd = Pandas.DataFrame(sort!(sceBhatt, :Class))

pd = Pandas.melt(pd, value_vars  = ["WT_UT__fm", "WT_2H__fm", "WT_8H__fm"], value_name = y, id_vars = ["GeneID","Class"] )
pd["Time after LPS (h)"] = [replace(replace(replace(split(ii, "_")[2], "UT" => 0), "2H" => 2), "8H" => 8) for ii in pd["variable"]]
pd["Genotype"] = [replace(split(ii, "_")[1], "RAD21" =>"Rad21KO") for ii in pd["variable"]]

Seaborn.boxplot(data = pd, y = y, x = "Class",hue = "Time after LPS (h)",showfliers = false,  palette = "Blues")
pretty_axes2()

title("WT")

hs = [  1, 1.1, 1.2,1, 1.1, 1.2, 1.1, 1.2, 1.3, 1, 1.1, 1.2, 1, 1.1, 1.2, 1.1, 1.2, 1.3, 1, 1.1, 1.2].*4 .*1.5
ps = adjust([ii for ii in p1as], Bonferroni())
u = 0.05 .*1.5

global iii = 0
for pp in 1:6
    global iii +=1 
    plt.plot([pp-1.25, pp-1.05], [hs[iii], hs[iii]], lw = 0.75, c= "black")
    plt.annotate("P = "*string(round(ps[iii], sigdigits = 2)),xy = [Statistics.mean([pp-1.25, pp-1.05]), hs[iii]+u], va = "center", ha = "center")

    global iii+=1
    plt.plot([pp-0.95, pp-0.75], [hs[iii], hs[iii]], lw = 0.75, c= "black")
    annotate("P = "*string(round(ps[iii], sigdigits = 2)),xy = [Statistics.mean([pp-0.95, pp-0.75]), hs[iii]+u], va = "center", ha = "center")

    global iii+=1
    plt.plot([pp-1.25, pp-0.75], [hs[iii], hs[iii]], lw = 0.75, c= "black")
    plt.annotate("P = "*string(round(ps[iii], sigdigits = 2)),xy = [Statistics.mean([pp-1.25, pp-0.75]), hs[iii]+u], va = "center", ha = "center")


end


squareplot()

legend_removal()

subplot(1,2,2)
y = "moment burst size"
pd = Pandas.DataFrame(sort!(sceBhatt, :Class))

pd = Pandas.melt(pd, value_vars  = ["WT_UT__bm", "WT_2H__bm", "WT_8H__bm"], value_name = y, id_vars = ["GeneID","Class"] )
pd["Time after LPS (h)"] = [replace(replace(replace(split(ii, "_")[2], "UT" => 0), "2H" => 2), "8H" => 8) for ii in pd["variable"]]
pd["Genotype"] = [replace(split(ii, "_")[1], "RAD21" =>"Rad21KO") for ii in pd["variable"]]

Seaborn.boxplot(data = pd, y = y, x = "Class",hue = "Time after LPS (h)",showfliers = false,  palette = "Blues")
pretty_axes2()



hs = [ 8.5, 8.75, 9, 8.5, 8.75, 9, 9.1, 9.35, 9.65, 9, 9.25, 9.5, 8.5, 8.75, 9, 9.1, 9.35, 9.65, 9, 9.25, 9.5].*500
ps = adjust([ii for ii in p1s], Bonferroni())
u = 0.09 * 500

global ii = 0

for pp in 1:6
    global ii+=1
    plt.plot([pp-1.25, pp-1.05], [hs[ii], hs[ii]], lw = 0.75, c= "black")
    plt.annotate("P = "*string(round(ps[ii], sigdigits = 2)),xy = [Statistics.mean([pp-1.25, pp-1.05]), hs[ii]+u], va = "center", ha = "center")

    global ii+=1
    plt.plot([pp-0.95, pp-0.75], [hs[ii], hs[ii]], lw = 0.75, c= "black")
    plt.annotate("P = "*string(round(ps[ii], sigdigits = 2)),xy = [Statistics.mean([pp-0.95, pp-0.75]), hs[ii]+u], va = "center", ha = "center")

    global ii+=1
    plt.plot([pp-1.25, pp-0.75], [hs[ii], hs[ii]], lw = 0.75, c= "black")
    plt.annotate("P = "*string(round(ps[ii], sigdigits = 2)),xy = [Statistics.mean([pp-1.25, pp-0.75]), hs[ii]+u], va = "center", ha = "center")

    
    

end

squareplot()
legend_out_of_plot()

plt.tight_layout()

title("WT")

savefigwithtext("BhattClass_WT_6class"*"percent"*string(alpha*100)*"__moments.svg")