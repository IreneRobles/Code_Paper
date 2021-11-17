

function Il12b_maturefig(gen)
    reps = CSV.read("../CompleteSets/GeneData/Il12b_intron_original.csv", DataFrames.DataFrame);
    reps = reps[reps[!,:Genotype].==gen, :]

    pd = Pandas.DataFrame(reps)
    Seaborn.boxplot(data = pd, x= "Timepoint", y = "N_exon",palette = "Blues", showfliers = false, )
     pvals = CSV.read("TukeyHSD_Il12b_mature_"*gen*".csv", DataFrames.DataFrame)
    #Seaborn.stripplot(data = pd, x= "Timepoint", y = "N_exon", hue = "Rep",palette = "Greys", size = 1.5,jitter = 0.45, zorder = 0)
    pretty_axes2()
    ylim(0, 70)

    h = 10
    u = 2

    xy = [0, 1]
    plt.plot(xy, [h, h], c = "black")
    p =   round(pvals[pvals[!,:Column1] .== gen*"_60-"*gen*"_0", "p adj"][1], sigdigits = 2)
  
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    h = 30
    xy = [1, 2]
    plt.plot(xy, [h, h], c = "black")
    p =   round(pvals[pvals[!,:Column1] .== ""*gen*"_90-"*gen*"_60", "p adj"][1], sigdigits = 2)
  
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    h = 60
    xy = [2, 3]
    plt.plot(xy, [h, h], c = "black")
    p =   round(pvals[pvals[!,:Column1] .== gen*"_90-"*gen*"_120", "p adj"][1], sigdigits = 2)
  
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    ylabel("Il12b  mature \n mRNA counts")
    xlabel("Time after LPS (min)")
    
    PrettyPlotting.line075black()
end


function Il12b_BF(gen; pal = "Blues")
    
    gen = "WT"
    s = MultipleTesting.adjust([do_mantelhaen(summary1, gen*"_0", gen*"_60")[1],
        do_mantelhaen(summary1, gen*"_60", gen*"_90")[1],
        do_mantelhaen(summary1, gen*"_90", gen*"_120")[1],
        ], Bonferroni())

    pd = Pandas.DataFrame(summary1[summary1[!,:Genotype].==gen, :])
    Seaborn.boxplot(data = pd, x= "Timepoint", y = "BurstFraction", palette = pal, showfliers = false, )
   # Seaborn.stripplot(data = pd, x= "Timepoint", y = "BurstFraction", hue = "Rep",palette = "Greys", size = 5,jitter = 0.1, zorder = 1)
    pretty_axes2()

    h = 0.05
    u = 0.005

    xy = [0, 1]
    plt.plot(xy, [h, h], c = "black")
    p =   round(s[1], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")
    h = 0.1
    xy = [1, 2]
    plt.plot(xy, [h, h], c = "black")
    p =   round(s[2], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")
    h = 0.12
    xy = [2, 3]
    plt.plot(xy, [h, h], c = "black")
    p =   round(s[3], sigdigits = 2)
    plt.annotate("$p", xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")
    ylabel("Il12b \n Burst Frequency")
    xlabel("Time after LPS (min)")
     PrettyPlotting.line075black()
end
