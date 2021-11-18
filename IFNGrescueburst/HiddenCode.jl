ENV["Code"] = "../../Code"
for folder in readdir(ENV["Code"]); push!(LOAD_PATH, normpath(ENV["Code"], folder));end

using CSV, RCall
using NoLongerProblems_FileHandling
using NoLongerProblems
using DataFrames
using FQfiles
using HypothesisTests
using MultipleTesting
using Seaborn
using Statistics
using PyPlot
using Random
using PrettyPlotting
using DataFrames



ENV["COLUMNS"] = 100000


function make_sumaries(tb, name_; int = "r2", limit = 1)
tb = tb[!, filter(c -> count(ismissing, tb[:,c])/size(tb,1) < 0.1, names(tb))]
tb[!,:N_intron] = tb[!, :N_thres_Nuc]
tb[!,:N_exon] = tb[!, :N_thres_Total]
REPS = [tb[tb[!,:Rep].== ii, :] for ii in unique(tb[!,:Rep])]
summary_= mean_burst_size_and_burst_fraction(REPS..., int = int, limit = limit)
CSV.write("summary_"*name_*".csv", summary_)
summary_
end

function BF_fig_90(summary1; pal = ["darkgray","red"], hs = [0.03, 0.05, 0.12, 0.05], u = maximum(hs)/20, tim = 90)

    s = MultipleTesting.adjust([do_mantelhaen(summary1, "WT_0", "Rad21KO_0")[1],
        do_mantelhaen(summary1, "WT_i0", "Rad21KO_i0")[1],
        do_mantelhaen(summary1, "WT_"*string(tim), "Rad21KO_"*string(tim))[1],
        do_mantelhaen(summary1, "WT_i"*string(tim), "Rad21KO_i"*string(tim))[1]
        ], Bonferroni())

    pd = Pandas.DataFrame(summary1)
    Seaborn.boxplot(data = pd, hue = "Genotype", x= "Timepoint", y = "BurstFraction", palette = pal, showfliers = false, )
   # Seaborn.stripplot(data = pd, x= "Timepoint", y = "BurstFraction", hue = "Rep",palette = "Greys", size = 5,jitter = 0.1, zorder = 1)
    pretty_axes2()

    h = hs[1]
    
    xy = [-0.25, 0.25]
    plt.plot(xy, [h, h], c = "black")
    p =   round(s[1], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")
    
    h = hs[2]
    xy = xy .+ 1
    plt.plot(xy, [h, h], c = "black")
    p =   round(s[2], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")
    
    h = hs[3]
    xy = xy .+ 1
    plt.plot(xy, [h, h], c = "black")
    p =   round(s[3], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")
    
     h = hs[4]
    xy = xy .+ 1
    plt.plot(xy, [h, h], c = "black")
    p =   round(s[4], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")
    
    
        ylim(0, maximum(hs)*1.2)

    ylabel("Burst Frequency")
    xlabel("Time after LPS (min)")
     PrettyPlotting.line075black()
    PrettyPlotting.squareplot()
    PrettyPlotting.legend_removal()
    
end


function do_bs_tests(df, name; limit = 1)
REPS = [df[df[!,:Rep].==ii, :] for ii in unique(df[!,:Rep])]
locusd = locusdata2(REPS...)

locusd = locusd[locusd[!,:BurstSize].> limit, :]
fname = "TukeyHSD_"*name*"_BurstSize.csv"
R"""
tb <- $locusd
tb$Timepoint <- as.factor(tb$Timepoint)
tb$Genotype <- as.factor(tb$Genotype)
tb$Sample <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb
s <- tb
aov.s = aov(BurstSize ~ Sample + Rep  ,data=s)  #do the analysis of variance
test <- TukeyHSD(aov.s)
test$Sample
write.csv(test$Sample,$fname)
test$Sample
"""
end

function Fig_BS(df, genname; pal = ["darkgray", "red"], n = 5, limit = 1, tim = 90)
    tim = string(tim)
    REPS = [df[df[!,:Rep].==ii, :] for ii in unique(df[!,:Rep])]
    locusd = locusdata2(REPS...)

    locusd = locusd[locusd[!,:BurstSize].> limit, :]
    
    pvals = CSV.read("TukeyHSD_"*genname*"_BurstSize.csv", DataFrames.DataFrame)


    pd = Pandas.DataFrame(locusd)
    Seaborn.boxplot(data = pd,hue = "Genotype", x= "Timepoint", y = "BurstSize", palette = pal, showfliers = false, )
    #Seaborn.stripplot(data = pd, x= "Timepoint", y = "BurstSize", hue = "Rep",palette = "Greys", size = 1.5,jitter = 0.45, zorder = 0)
    pretty_axes2()
    ylim(0, 22*n)

    h = 19.5*n
    u = 0.7*n

    xy = [-0.25, 0.25]
    plt.plot(xy, [h, h], c = "black")
    p =   round(pvals[.|(pvals[!,:Column1] .== "Rad21KO_0-WT_0", pvals[!,:Column1] .== "WT_0-Rad21KO_0"), "p adj"][1], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")
    
    
    h = 19.5*n*1.05
    xy =xy .+ 1
    plt.plot(xy, [h, h], c = "black")
    p =   round(pvals[.|(pvals[!,:Column1] .== "Rad21KO_i0-WT_i0", pvals[!,:Column1] .== "WT_i0-Rad21KO_i0"), "p adj"][1], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    
    h = 19.5*n*1
    xy =xy .+ 1
    plt.plot(xy, [h, h], c = "black")
    p =   round(pvals[.|(pvals[!,:Column1] .== "Rad21KO_"*tim*"-WT_"*tim, pvals[!,:Column1] .== "WT_"*tim*"-Rad21KO_"*tim), "p adj"][1], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    
    
    
    h = 19.5*n*1
    xy =xy .+ 1
    plt.plot(xy, [h, h], c = "black")
    p =   round(pvals[.|(pvals[!,:Column1] .== "Rad21KO_i"*tim*"-WT_i"*tim, pvals[!,:Column1] .== "WT_i"*tim*"-Rad21KO_i"*tim), "p adj"][1], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    
    
    ylabel("$name \n Burst Size")
    xlabel("Time after LPS (min)")
    
        PrettyPlotting.line075black()

end

