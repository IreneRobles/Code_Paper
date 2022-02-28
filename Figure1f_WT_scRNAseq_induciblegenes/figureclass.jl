figure(figsize = (20, 10))


subplot(1,2,1)
y = "Fraction of cells \n with transcripts detected"
pd = Pandas.DataFrame(sort!(sceBhatt, :Class))

pd = Pandas.melt(pd, value_vars  = ["WT_UT__alpha", "WT_2H__alpha", "WT_8H__alpha"], value_name = y, id_vars = ["GeneID","Class"] )
pd["Time after LPS (h)"] = [replace(replace(replace(split(ii, "_")[2], "UT" => 0), "2H" => 2), "8H" => 8) for ii in pd["variable"]]
pd["Genotype"] = [replace(split(ii, "_")[1], "RAD21" =>"Rad21KO") for ii in pd["variable"]]

Seaborn.boxplot(data = pd, y = y, x = "Class",hue = "Time after LPS (h)",showfliers = false,  palette = "Blues")
pretty_axes2()

title("WT")

hs = [  1, 1.1, 1.2,1, 1.1, 1.2, 1.1, 1.2, 1.3, 1, 1.1, 1.2, 1, 1.1, 1.2, 1.1, 1.2, 1.3, 1, 1.1, 1.2]
ps = adjust([ii for ii in p1as], Bonferroni())
global iii = 0
for pp in 1:6
    global iii +=1 
    plot([pp-1.25, pp-1.05], [hs[iii], hs[iii]], lw = 0.75, c= "black")
    annotate("P = "*string(round(ps[iii], sigdigits = 2)),xy = [Statistics.mean([pp-1.25, pp-1.05]), hs[iii]+0.05], va = "center", ha = "center")

    global iii+=1
    plot([pp-0.95, pp-0.75], [hs[iii], hs[iii]], lw = 0.75, c= "black")
    annotate("P = "*string(round(ps[iii], sigdigits = 2)),xy = [Statistics.mean([pp-0.95, pp-0.75]), hs[iii]+0.05], va = "center", ha = "center")

    global iii+=1
    plot([pp-1.25, pp-0.75], [hs[iii], hs[iii]], lw = 0.75, c= "black")
    annotate("P = "*string(round(ps[iii], sigdigits = 2)),xy = [Statistics.mean([pp-1.25, pp-0.75]), hs[iii]+0.05], va = "center", ha = "center")


end

ylim(0, 1.5)
squareplot()

legend_removal()

subplot(1,2,2)
y = "ln(CPM + 1) \n cells with transcripts detected"
pd = Pandas.DataFrame(sort!(sceBhatt, :Class))

pd = Pandas.melt(pd, value_vars  = ["WT_UT__mu", "WT_2H__mu", "WT_8H__mu"], value_name = y, id_vars = ["GeneID","Class"] )
pd["Time after LPS (h)"] = [replace(replace(replace(split(ii, "_")[2], "UT" => 0), "2H" => 2), "8H" => 8) for ii in pd["variable"]]
pd["Genotype"] = [replace(split(ii, "_")[1], "RAD21" =>"Rad21KO") for ii in pd["variable"]]

Seaborn.boxplot(data = pd, y = y, x = "Class",hue = "Time after LPS (h)",showfliers = false,  palette = "Blues")
pretty_axes2()



hs = [ 8.5, 8.75, 9, 8.5, 8.75, 9, 9.1, 9.35, 9.65, 9, 9.25, 9.5, 8.5, 8.75, 9, 9.1, 9.35, 9.65, 9, 9.25, 9.5]
ps = adjust([ii for ii in p1s], Bonferroni())

global ii = 0

for pp in 1:6
    global ii+=1
    plot([pp-1.25, pp-1.05], [hs[ii], hs[ii]], lw = 0.75, c= "black")
    annotate("P = "*string(round(ps[ii], sigdigits = 2)),xy = [Statistics.mean([pp-1.25, pp-1.05]), hs[ii]+0.09], va = "center", ha = "center")

    global ii+=1
    plot([pp-0.95, pp-0.75], [hs[ii], hs[ii]], lw = 0.75, c= "black")
    annotate("P = "*string(round(ps[ii], sigdigits = 2)),xy = [Statistics.mean([pp-0.95, pp-0.75]), hs[ii]+0.09], va = "center", ha = "center")

    global ii+=1
    plot([pp-1.25, pp-0.75], [hs[ii], hs[ii]], lw = 0.75, c= "black")
    annotate("P = "*string(round(ps[ii], sigdigits = 2)),xy = [Statistics.mean([pp-1.25, pp-0.75]), hs[ii]+0.09], va = "center", ha = "center")

    
    

end

squareplot()
legend_out_of_plot()

plt.tight_layout()

title("WT")

savefigwithtext("BhattClass_WT_6class"*"percent"*string(alpha*100)*".svg")