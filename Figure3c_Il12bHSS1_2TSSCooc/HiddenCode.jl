ENV["Code"] = "../../Code"
for folder in readdir(ENV["Code"]); push!(LOAD_PATH, normpath(ENV["Code"], folder));end

using CSV
using NoLongerProblems_FileHandling
using NoLongerProblems
using NoLongerProblems_Pandas
using DataFrames
using FQfiles
using HypothesisTests
using MultipleTesting
using Seaborn
using Random
import Pandas
using Statistics, PyPlot
using PrettyPlotting
using DataFrames
include("../TSS_quantification/StandardCode.jl")
using RCall
using Distributions


function linked_data(genefolder,nascent,ehn,suff) 
    tb1 = DataFrame(CSV.File(normpath(genefolder, nascent*"__"*ehn*"__linked"*suff*".csv")))
    sort!(tb1, :Image_Cell)
    suff2 = if endswith(suff, "noz") "200nm_noz" else "200nm" end
    tb2 = DataFrame(CSV.File(normpath(genefolder, nascent*"__"*ehn*"__linked"*suff2*".csv")))
    sort!(tb2, :Image_Cell)
    # filters are only important when there is 1 and 1 burst, otherwise they should be linked
    bool1 =  tb2[!,nascent*"_N"].== 1
    bool2 =  tb1[!,ehn*"_N"].== 1
    tt1 = tb1[bool1.*bool2, :]
    tt2 = tb2[.! (bool1.*bool2), :]
    tb = vcat(tt1, tt2)
    
    for ii in  columns_containing(tb, nascent)
        rename!(tb, ii => replace(ii, nascent => :Gene))
    end
     for ii in  columns_containing(tb, ehn)
        rename!(tb, ii => replace(ii, ehn => :Enh))
    end
   tb
end

function get_pairs_with_distances(tb; limit_gene = 1, limit_enh =1)
    locus1 = tb[tb[!,:locus1_Gene_Enh].!="NA", :]
    cols = ["Image_Cell", "AREA_cell", "AREA_nuc", "Genotype", "Timepoint", "Rep", "Sample", "Cell",
        "locus2_Gene", "locus2_Enh", "locus2_Gene_Enh", "Gene_N", "Enh_N", "Image", 
        "locus2_Gene_size", "locus2_Enh_size"]
    locus2 = tb[tb[!,:locus2_Gene_Enh].!="NA", cols]
    
     for ii in  columns_containing(locus2, "locus2")
        rename!(locus2, ii => replace(ii, "locus2" => "locus1"))
    end
    tb = join_in_all_common_columns(locus1, locus2)
    tb = tb[tb[!,"locus1_Gene_size"].> limit_gene, :]
     tb = tb[tb[!,"locus1_Enh_size"].> limit_enh, :]
    return tb

end





function show_distances(genefolder,nascent,ehn,suff; limit = 4, limit_gene = 1, limit_enh =1) 
    tb = linked_data(genefolder,nascent,ehn,suff) 
    tb = tb[tb[!,:Timepoint].>0, :]
    if genefolder == "Prdm1"
        tb = tb[tb[!,:Timepoint].==60, :]
    end
    pairs = get_pairs_with_distances(tb)
        pairs = pairs[parse.(Float64,pairs[!,:locus1_Gene_Enh]).< limit, :]
    sort!(pairs, [:Genotype, :Timepoint], rev = true)

    pd = Pandas.DataFrame(pairs)
    Seaborn.boxplot(data = pd, y = "locus1_Gene_Enh", x = "Genotype", showfliers = false, palette = ["darkgray", "red"])
    pretty_axes2()
    title(nascent)
    legend_removal()
    ylabel("Distance Gene-Enhancer (um)")
    line075black()
    ylim(0, limit +2)
    xticks(rotation = 0)
    squareplot()
    
    
t = R"""


tb <- $pairs

tb$Genotype <- as.factor(tb$Genotype)
tb$Timepoint <- as.factor(tb$Timepoint)
tb$Rep <- as.factor(tb$Rep);

s <- tb
aov.s = aov(locus1_Gene_Enh ~ Genotype + Rep  ,data=s)  #do the analysis of variance
test <- TukeyHSD(aov.s)
test$Genotype
    """
       plot([0,1], [limit+limit/10, limit+limit/10], c = "black", lw = 0.75)
    annotate(string(round(t[4], sigdigits = 2 )), xy = [0.5, limit+limit/6],va = "center", ha = "center")
    
    return tb, pairs
    
end

function coocurrences(genefolder, nascent, enh; suff = "2nm_noz")
    tb = DataFrame(CSV.File(normpath(genefolder, nascent*"__"*enh*"__coocurrences"*suff*".csv")))
end


function OneTwoTSS(genefolder, target)
    tb = DataFrame(CSV.File(normpath(genefolder, target*".csv")))
    tb[!,:Sample_Rep] = tb[!,:Sample] .* "_" .* string.(tb[!,:Rep])
    sam_rep = unique(tb[!,:Sample_Rep])
    new_df = DataFrames.DataFrame(
    Sample_Rep = sam_rep,
    Genotype = [split(ii, "_")[1] for ii in sam_rep],
    Timepoint = [split(ii, "_")[2] for ii in sam_rep],
    Rep = [split(ii, "_")[3] for ii in sam_rep],
    N_cells = [sum(tb[!,:Sample_Rep].==ii) for ii in sam_rep],
    TSS_1 = [sum(tb[tb[!,:Sample_Rep].==ii, :N_TSS2].==1) for ii in sam_rep],
    TSS_2 = [sum(tb[tb[!,:Sample_Rep].==ii, :N_TSS2].==2) for ii in sam_rep]
    )
    new_df[!,:BF] = (new_df[!,:TSS_1] + 2 .* new_df[!,:TSS_2])./ (2 .*new_df[!,:N_cells])
    new_df[!,:TSS_1_F] = new_df[!,:TSS_1] ./new_df[!,:N_cells]
    new_df[!,:TSS_1_Frandom] = (2 .* new_df[!,:BF]) .- 2 .* (new_df[!,:BF].^2)
    new_df[!,:TSS_1_random] = round.(new_df[!,:TSS_1_Frandom].*new_df[!,:N_cells], digits = 0)
    new_df[!,:TSS_2_F] = new_df[!,:TSS_2] ./new_df[!,:N_cells]
    new_df[!,:TSS_2_Frandom] =  new_df[!,:BF].^2
    new_df[!,:TSS_2_random] = round.(new_df[!,:TSS_2_Frandom].*new_df[!,:N_cells], digits = 0)
    CSV.write( normpath(genefolder, target* "_OneTwoTSS.csv"),new_df)
    new_df
end

function OneTwoTSS_plot(genefolder, target, genotype; h1 = 0.01, h2 = 0.05)
    onetwo = OneTwoTSS(genefolder, target)
    onetwo = onetwo[onetwo[!,:Genotype].==genotype, :]
      if genefolder == "Prdm1"
        bool1 = onetwo[!,:Timepoint].=="60"
        bool2 = onetwo[!,:Timepoint].=="0"
         onetwo = onetwo[bool1.|bool2, :]
    end
    pd = Pandas.DataFrame(onetwo)
    
    pd1 = Pandas.melt(pd, value_vars = ["TSS_2", "TSS_2_random"], id_vars = [:Rep, :Timepoint, :Genotype, :N_cells], 
        value_name = "Count", var_name = "Observed_Random")
    pd1[:yes_no] = "yes"
    
    pd2 = Pandas.melt(pd, value_vars = ["TSS_2", "TSS_2_random"], id_vars = [:Rep, :Timepoint, :Genotype, :N_cells], 
        value_name = "Count", var_name = "Observed_Random")
    pd2[:Count] = pd2[:N_cells] .- pd2[:Count]
    pd2[:yes_no] = "no"

    
    onetwo_melt =vcat( DataFrames.DataFrame(pd1), DataFrames.DataFrame(pd2))
    onetwo_melt[!,:Observed_Random] = [if ii == "TSS_2" "Observed" else "Expected" end for ii in onetwo_melt[!,:Observed_Random]]
    onetwo_melt
    pvals = []
    
    for t in  sort!(unique(onetwo_melt[!,:Timepoint]))
         time = onetwo_melt[onetwo_melt[!,:Timepoint].==t, :]
        test = R"""
        library("dplyr")
        tb = $time
        tb$Count = as.integer(tb$Count) 

        Data <- mutate(tb,
           Observed_Random = factor(Observed_Random, levels=unique(Observed_Random)),
           yes_no = factor(yes_no, levels=unique(yes_no)),
           Rep = factor(Rep, levels=unique(Rep))
           )

# Last variable is the strata (the variable that is not check for assotiation)
Data.xtabs <- xtabs(Count ~ Observed_Random + yes_no + Rep, 
                       data=Data)

ftable(Data.xtabs)  

mantelhaen.test(Data.xtabs)
"""
         push!(pvals, test[3][1])
    end
    
    pd3 = Pandas.melt( Pandas.DataFrame(onetwo), value_vars = ["TSS_2_F", "TSS_2_Frandom"], id_vars = [:Rep, :Timepoint, :Genotype, :N_cells], 
        value_name = "Fraction of Cells", var_name = "2 TSS detected") 
    pd3["2 TSS detected"] = [if ii == "TSS_2_F" "Observed" else "Expected" end for ii in pd3["2 TSS detected"]]

    Seaborn.boxplot(data = pd3, y = "Fraction of Cells", x = "Timepoint", hue = "2 TSS detected", showfliers = false, palette = ["gray", "lightgray"])
    
    

    
    plot([-0.25,0.25], [h1, h1], c = "black", lw = 0.75)
    annotate(string(round(pvals[1], sigdigits = 3 )), xy = [0, h1 + h2/10],va = "center", ha = "center")
    
    plot([0.75,1.25], [h2, h2], c = "black", lw = 0.75)
    annotate(string(round(pvals[2], sigdigits = 3 )), xy = [1, h2 + h2/10],va = "center", ha = "center")
    
    ylim(0, h2*1.4)
    pretty_axes2()
    squareplot()
    line075black()
    ylabel("Fraction of cells with 2 TSSs")
title(target*" ($genotype)")
    legend_removal()
end