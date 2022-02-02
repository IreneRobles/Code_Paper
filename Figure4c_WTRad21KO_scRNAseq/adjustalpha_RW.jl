
sce = SingleCellExperiment.fit_mu_std_alpha(sce, splitdataby = :Sample, assay = "lnCPMplus1")
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
subsubsce = SingleCellExperiment.Shalek2014_module_score(collect(bhattgenes[!,"GeneID"]), subsubsce,fitparameter = "mu", modulescore_name = :BhattGenesScore, untreated_pattern = "UT",comparedtothissample = "WT", assay = "CPM")

sceBhatt = innerjoin(subsubsce.rowData, bhattgenes,on = :GeneID)
sceBhatt[!,:Class] = [replace(replace(ii, "A1"=>"A1+2"), "A2"=>"A1+2") for ii in  sceBhatt[!,:Class]]

println(string("Percent ", alpha*100, "%"))
println(string("Total genes considered ", nrow(sceBhatt)))

#println("ln CPM plus 1 in expressing cells")
#println("WT UT vs Rad21KO UT")

t = HypothesisTests.SignedRankTest([ii for ii in sceBhatt[!,"WT_UT__mu"]], [ii for ii in sceBhatt[!,"RAD21_UT__mu"]])
#println(t)
p1 = pvalue(t)


#println("Fraction Expressing cells")
#println("WT UT vs Rad21KO UT")

t = HypothesisTests.SignedRankTest([ii for ii in sceBhatt[!,"WT_UT__alpha"]], [ii for ii in sceBhatt[!,"RAD21_UT__alpha"]])
#println(t)
p1a = pvalue(t)

#println("ln CPM plus 1 in expressing cells")
#println("WT 2H vs Rad21KO 2H")

t = HypothesisTests.SignedRankTest([ii for ii in sceBhatt[!,"WT_2H__mu"]], [ii for ii in sceBhatt[!,"RAD21_2H__mu"]])
#println(t)
p2 = pvalue(t)

#println("Fraction Expressing cells")
#println("WT 2H vs Rad21KO 2H")

t = HypothesisTests.SignedRankTest([ii for ii in sceBhatt[!,"WT_2H__alpha"]], [ii for ii in sceBhatt[!,"RAD21_2H__alpha"]])
#println(t)
p2a = pvalue(t)

#println("ln CPM plus 1 in expressing cells")
#println("WT 8H vs Rad21KO 8H")
t = HypothesisTests.SignedRankTest([ii for ii in sceBhatt[!,"WT_8H__mu"]], [ii for ii in sceBhatt[!,"RAD21_8H__mu"]])
#println(t)
p3 = pvalue(t)
#println("Fraction Expressing cells")

#println("WT 8H vs Rad21KO 8H")
t = HypothesisTests.SignedRankTest([ii for ii in sceBhatt[!,"WT_8H__alpha"]], [ii for ii in sceBhatt[!,"RAD21_8H__alpha"]])
#println(t)
p3a = pvalue(t)



figure(figsize = (7, 6))


subplot(1,2,2)
y = "ln(CPM + 1) \n cells with transcripts detected"
pd = Pandas.DataFrame(sort!(sceBhatt, :Class))

pd = Pandas.melt(pd, value_vars  = ["WT_UT__mu", "RAD21_UT__mu" ,"WT_2H__mu", "RAD21_2H__mu", "WT_8H__mu", "RAD21_8H__mu"], value_name = y, id_vars = ["GeneID","Class"] )
pd["Time after LPS (h)"] = [replace(replace(replace(split(ii, "_")[2], "UT" => 0), "2H" => 2), "8H" => 8) for ii in pd["variable"]]
pd["Genotype"] = [replace(split(ii, "_")[1], "RAD21" =>"Rad21KO") for ii in pd["variable"]]

Seaborn.boxplot(data = pd, y = y,x = "Time after LPS (h)", hue = "Genotype",showfliers = false,  palette = ["darkgray", "red"])
pretty_axes2()

hs = [8.2, 8.7, 9.0]
ps = adjust([p1,p2,p3], Bonferroni())

for ii in 1:length(hs)
plt.plot([-0.25+ii-1, 0.25+ii-1], [hs[ii], hs[ii]], lw = 0.75, c= "black")
annotate("P = "*string(round(ps[ii], sigdigits = 3)),xy = [ii-1, hs[ii]+0.1], va = "center", ha = "center")
end
squareplot()

ax = gca()
for line in ax.get_lines()
    line.set_color("black")
end
subplot(1,2,1)



y = "Fraction of cells \n with transcripts detected"
pd = Pandas.DataFrame(sort!(sceBhatt, :Class))

pd = Pandas.melt(pd, value_vars  = ["WT_UT__alpha", "RAD21_UT__alpha" ,"WT_2H__alpha", "RAD21_2H__alpha", "WT_8H__alpha", "RAD21_8H__alpha"], value_name = y, id_vars = ["GeneID","Class"] )
pd["Time after LPS (h)"] = [replace(replace(replace(split(ii, "_")[2], "UT" => 0), "2H" => 2), "8H" => 8) for ii in pd["variable"]]
pd["Genotype"] = [replace(split(ii, "_")[1], "RAD21" =>"Rad21KO") for ii in pd["variable"]]

Seaborn.boxplot(data = pd, y = y,x = "Time after LPS (h)", hue = "Genotype",showfliers = false, palette = ["darkgray", "red"])
pretty_axes2()

hs = [0.6, 1.05, 1.05]
ps = adjust([p1a,p2a,p3a], Bonferroni())

for ii in 1:length(hs)
plt.plot([-0.25+ii-1, 0.25+ii-1], [hs[ii], hs[ii]], lw = 0.75, c= "black")
annotate("P = "*string(round(ps[ii], sigdigits = 2)),xy = [ii-1, hs[ii]+0.05], va = "center", ha = "center")
end

legend_removal()

ylim(-0.05, 1.1)
squareplot()
ax = gca()
for line in ax.get_lines()
    line.set_color("black")
end
plt.tight_layout()
savefigwithtext("scRNAseq_mu_alpha_bhattgenes_WTRad21KO_percent"*string(alpha*100)*".svg")

