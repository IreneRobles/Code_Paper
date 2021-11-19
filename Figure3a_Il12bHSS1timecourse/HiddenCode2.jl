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
using PyCall
using Statistics
using PyPlot
using PrettyPlotting
using DataFrames

function plotexons_violin(exp; y = "N_exon")
    pd = Pandas.DataFrame(exp);
    Seaborn.violinplot(data = pd, y = y, x = "Timepoint", hue = "Genotype", cut = 0)
    Seaborn.stripplot(data = pd, y = y, x = "Timepoint", hue = "Genotype", dodge = true, zorder = 0)
    PrettyPlotting.pretty_axes2()
end


function plotexons_box(exp;y = "N_exon")
    pd = Pandas.DataFrame(exp);
    Seaborn.boxplot(data = pd, y = y, x = "Timepoint", hue = "Genotype")
    Seaborn.stripplot(data = pd, y = y, x = "Timepoint", hue = "Genotype", dodge = true, zorder = 0)
    PrettyPlotting.pretty_axes2()
end

function cells_per_sample(df; col = "Sample", newcol= "N_Cells")
    sams = unique(df[!,col])
    ncells = [ count(x -> x == sam, df[!,col]) for sam in sams]
    new_df = DataFrames.DataFrame(:Sample=> sams, Symbol(newcol)=> ncells)
end


function fraction_positive(df, column, limit; col = "Sample", newcol= "Positive")
    sams = unique(df[!,col])
    ncells = [ count(x -> x == sam, df[!,col]) for sam in sams]
    sp_df = [ df[df[!,col].==sam, :] for sam in sams]
    pos_cells = [ sum(ii[!,column].>limit) for ii in sp_df]
    new_df = DataFrames.DataFrame(:Sample=> sams, Symbol("N_Cells")=> ncells, Symbol(newcol*"$limit")=> pos_cells,
    Symbol(newcol*"$limit"*"_Fraction")=> pos_cells./ncells)
end



function fraction_positive_plot(df, column, SAMPLE; col = "Sample", color = "black", max = 20, addtosam = "")
    subdf = df[df[!,col].==SAMPLE, :]
    plot(0:max, [Statistics.mean(subdf[!,column].>ii) for ii in 0:max], c = color, label = SAMPLE*addtosam)
    ylabel("Fraction of cells with \n more than")
    xlabel("counts")
end


function fraction_positive_plot_wt(exp; max = 20, y = "N_exon")

s = "WT_0"
fraction_positive_plot(exp, y, s; col = "Sample", color = "#87CEEB", max = max, addtosam = "")
s = "WT_90"
fraction_positive_plot(exp, y, s; col = "Sample", color = "#1B3F8B", max = max, addtosam = "")
pretty_axes2()

end


function fraction_positive_plot_rad21(exp; max = 20, y = "N_exon")

s = "Rad21KO_0"
fraction_positive_plot(exp, y, s; col = "Sample", color = "#ffffb2", max = max, addtosam = "")

s = "Rad21KO_90"

    fraction_positive_plot(exp, y, s; col = "Sample", color = "#e31a1c", max = max, addtosam = "")
pretty_axes2()

end

function sample_from_genotype_timepoint(df)
    df[!,:Sample] = df[!,:Genotype] .*"_".*string.(df[!,:Timepoint])
    df
end


function calculatesumary(expe; r = 2, limit = 1)
    sams = unique(expe[!, :Sample])
    rep = unique(expe[!, :Rep])[1]
    sp_df = split_by(expe, :Sample)
    N_Cells = [nrow(sp_df[sam]) for sam in sams]
    N_TSS1 = [sum(sp_df[sam][!,Symbol("TSS1_r"*string(r))].> limit) for sam in sams]
    N_TSS2 = [sum(sp_df[sam][!,Symbol("TSS2_r"*string(r))].> limit) for sam in sams]
    new_df = DataFrames.DataFrame(
        :Sample =>sams,
        :Rep =>[rep for ii in sams],
        :Genotype =>[split(s, "_")[1] for s in sams],
        :Timepoint =>[split(s, "_")[2] for s in sams],
        :N_cells => N_Cells,
        :N_TSS1=>N_TSS1,
        :N_TSS2=>N_TSS2,
        :N_TSS=>N_TSS1.+N_TSS2,
        :BF =>(N_TSS1.+N_TSS2)./2N_Cells
    )
end


function calculatesumary(expe...; r = 2, limit = 1)
   join_in_all_common_columns([calculatesumary(ex; r = r, limit = limit)  for ex in expe]...)
end

function get_data(genefolder, genename, expnumber; orgainserfolder =  "../OrganiseData")
    expnumber = string(expnumber)
    folder = orgainserfolder * "/"*genefolder*"/" 
     exp = CSV.read(folder * genename*"_exp"*expnumber*"_FQ.csv", DataFrames.DataFrame)
    t = CSV.read(folder * genename*"_exp"*expnumber*"_TSS.csv", DataFrames.DataFrame)
    t[!,:Image_Cell] = t[!,:Image] .*"__".*string.(t[!,:Cell])
    exp = innerjoin(exp, t, on = :Image_Cell, makeunique=true)
    cp = CSV.read(folder * genename*"_exp"*expnumber*"_CP.csv", DataFrames.DataFrame)[!,["Genotype","Timepoint", "Rep","Image_Cell"]]
    expreturn = sample_from_genotype_timepoint(sort!(innerjoin(exp, cp, on = :Image_Cell, makeunique=true), :Genotype, rev = true))
    expreturn
end

function cells_per_replicate(bigdf; AnalysisName = AnalysisName)

    dfs = DataFrames.DataFrame()

    for ii in unique(bigdf[!,:Rep])
        sub = bigdf[bigdf[!,:Rep].== ii, :]
        coundf = cells_per_sample(sub, col = "Sample", newcol= "N_Cells_Rep"*string(ii))
        if ii == unique(bigdf[!,:Rep])[1]
            dfs = coundf
        else
            dfs = innerjoin(dfs, coundf, on = :Sample)
        end
    end

    CSV.write(AnalysisName*"/"*AnalysisName*"_ncells_replicate.csv")
    dfs
end


function plotexons_violin(exp; y = "N_exon")
    pd = Pandas.DataFrame(exp);
    Seaborn.violinplot(data = pd, y = y, x = "Timepoint", hue = "Genotype", cut = 0)
    Seaborn.stripplot(data = pd, y = y, x = "Timepoint", hue = "Genotype", dodge = true, zorder = 0)
    PrettyPlotting.pretty_axes2()
end


function plotexons_box(exp;y = "N_exon")
    pd = Pandas.DataFrame(exp);
    Seaborn.boxplot(data = pd, y = y, x = "Timepoint", hue = "Genotype")
    Seaborn.stripplot(data = pd, y = y, x = "Timepoint", hue = "Genotype", dodge = true, zorder = 0)
    PrettyPlotting.pretty_axes2()
end

function cells_per_sample(df; col = "Sample", newcol= "N_Cells")
    sams = unique(df[!,col])
    ncells = [ count(x -> x == sam, df[!,col]) for sam in sams]
    new_df = DataFrames.DataFrame(:Sample=> sams, Symbol(newcol)=> ncells)
end


function fraction_positive(df, column, limit; col = "Sample", newcol= "Positive", rep = "1")
    sams = unique(df[!,col])
    ncells = [ count(x -> x == sam, df[!,col]) for sam in sams]
    sp_df = [ df[df[!,col].==sam, :] for sam in sams]
    pos_cells = [ sum(ii[!,column].>limit) for ii in sp_df]
    new_df = DataFrames.DataFrame(Symbol(col)=> sams, Symbol("N_Cells")=> ncells, Symbol(newcol*"$limit")=> pos_cells,
    Symbol(newcol*"$limit"*"_Fraction")=> pos_cells./ncells)
    new_df[!,:Rep] = [rep for ii in 1:nrow(new_df)]
    new_df
end





function fraction_positive_plot(df, column, SAMPLE; col = "Sample", color = "black", max = 20, addtosam = "")
    subdf = df[df[!,col].==SAMPLE, :]
    plot(0:max, [Statistics.mean(subdf[!,column].>ii) for ii in 0:max], c = color, label = SAMPLE*addtosam)
    ylabel("Fraction of cells with \n more than")
    xlabel("counts")
end


function fraction_positive_plot_wt(exp; max = 20)

s = "WT_0"
fraction_positive_plot(exp, "N_exon", s; col = "Sample", color = "#87CEEB", max = max, addtosam = "")
s = "WT_120"
fraction_positive_plot(exp, "N_exon", s; col = "Sample", color = "#1B3F8B", max = max, addtosam = "")
pretty_axes2()

end


function fraction_positive_plot_rad21(exp; max = 20)

s = "Rad21KO_0"
fraction_positive_plot(exp, "N_exon", s; col = "Sample", color = "#ffffb2", max = max, addtosam = "")

s = "Rad21KO_120"
fraction_positive_plot(exp, "N_exon", s; col = "Sample", color = "#e31a1c", max = max, addtosam = "")
pretty_axes2()

end

function addrep(experiment, n)
    experiment[!,:Rep] = [n for ii in 1:nrow(experiment)]
    experiment
end


using RCall

function plot_wt_frequencies(GENE; activationtimepoint = 90)
    set = sort!(get_genedata(GENE), :Timepoint)
    if GENE == "Hprt" ||GENE == "Fh1"
        set = set[set[!, :Timepoint].!="0i", :]
        set = set[set[!, :Timepoint].!="8i", :]
        set[!, :Timepoint] = [parse(Int, ii)*60 for ii in set[!, :Timepoint]]
    end
    wt_data = set[set[!,:Genotype].=="WT", :]
    
    reps = split_by(wt_data, :Rep)
    wt = vcat([fraction_positive(reps[ii], "total_TS_Cell", 0; col = "Timepoint", newcol= "TSSMorethan", rep = ii) for ii in keys(reps)]...)
    wt_tss2 = vcat([fraction_positive(reps[ii], "total_TS_Cell", 1; col = "Timepoint", newcol= "TSSMorethan", rep = ii) for ii in keys(reps)]...)
    wt[!,"TSSMorethan1"] = wt_tss2[!,"TSSMorethan1"]
    wt[!,"TSSMorethan1_Fraction"] = wt_tss2[!,"TSSMorethan1_Fraction"]
    wt[!,"N_TSS2"] = wt[!,"TSSMorethan0"] .+ wt[!,"TSSMorethan1"]
    wt[!,"Burst Fraction"] = (wt[!,"TSSMorethan0"] .+ wt[!,"TSSMorethan1"]) ./ (2* wt[!,"N_Cells"])
    
    wt[!,:Genotype] = ["WT" for ii in 1:nrow(wt)]
    wt[!,:Sample] =   wt[!,:Genotype].*"_".*  string.(wt[!,:Timepoint])
    Seaborn.boxplot(data = Pandas.DataFrame(wt), x = "Timepoint", y = "Burst Fraction", palette = "Blues")
    Seaborn.stripplot(data = Pandas.DataFrame(wt), x = "Timepoint", y = "Burst Fraction", hue = "Rep", palette = "Set2")
    pretty_axes2()
    
    title(GENE)
    
    legend_removal()
    
    xlabel("Time after LPS (min)")
    
    return wt

end




function do_mantelhaen(df, comp1, comp2)
    s1 = df[df[!,:Sample] .== comp1, :]
    s2 = df[df[!,:Sample] .== comp2, :]
    tb_test = join_in_all_common_columns(s1, s2)
    tb1 = tb_test; tb1[!,:Burst] = ["Yes" for ii in 1:nrow(tb_test)]
    tb1[!,:Count] = tb1[!,:N_TSS2]
    tb1 = tb1[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb2 = tb_test; tb2[!,:Burst] = ["No" for ii in 1:nrow(tb_test)]
    tb2[!,:Count] = tb2[!,:N_Cells] * 2 - tb2[!,:N_TSS2]
    tb2 = tb2[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb = join_in_all_common_columns(tb1, tb2)

    
test = R"""
library("psych")
library("vcd")
library("DescTools")
library("rcompanion")
library("dplyr")

tb = $tb
tb$Count = as.integer(tb$Count)
    


    

Data <- mutate(tb,
           Sample = factor(Sample, levels=unique(Sample)),
           Burst = factor(Burst, levels=unique(Burst)),
           Rep = factor(Rep, levels=unique(Rep))
           )

# Last variable is the strata (the variable that is not check for assotiation)
Data.xtabs <- xtabs(Count ~ Burst + Sample + Rep, 
                       data=Data)

ftable(Data.xtabs)  

mantelhaen.test(Data.xtabs)
"""

return test[3][1],test[5][1]
end


function add_tests(tests, hi; u = 0.005)
    
pvals = adjust([t[1] for t in tests], Bonferroni())

for ii in 1:length(pvals)
    h = hi[ii]
    plot([ii-1, ii], [h, h], c = "Black")
    annotate(NoLongerProblems.transform_pvalue_in_stars(pvals[ii]), xy =[mean([ii-1, ii]), h+u], ha = "center", va = "center")
    
end
end


function add_tests2(pvals, hi; u = 0.005)
    
for ii in 1:length(pvals)
    h = hi[ii]
    plot([ii-1, ii], [h, h], c = "Black")
    annotate(NoLongerProblems.transform_pvalue_in_stars(pvals[ii]), xy =[mean([ii-1, ii]), h+u], ha = "center", va = "center")
    
end
end



function calculate_means(df, column; col = "Sample", newcol= "Mean", rep = "1")
    sams = unique(df[!,col])
    ncells = [ count(x -> x == sam, df[!,col]) for sam in sams]
    sp_df = [ df[df[!,col].==sam, :] for sam in sams]
    pos_cells = [ mean(ii[!,column]) for ii in sp_df]
     new_df = DataFrames.DataFrame(Symbol(col)=> sams, Symbol("N")=> ncells,
    Symbol(newcol)=> pos_cells)
    new_df[!,:Rep] = [rep for ii in 1:nrow(new_df)]
    new_df
end

function calculate_medians(df, column; col = "Sample", newcol= "Median", rep = "1")
    sams = unique(df[!,col])
    ncells = [ count(x -> x == sam, df[!,col]) for sam in sams]
    sp_df = [ df[df[!,col].==sam, :] for sam in sams]
    pos_cells = [ mean(ii[!,column]) for ii in sp_df]
     new_df = DataFrames.DataFrame(Symbol(col)=> sams, Symbol("N")=> ncells,
    Symbol(newcol)=> pos_cells)
    new_df[!,:Rep] = [rep for ii in 1:nrow(new_df)]
    new_df
end


function tss_data(tss; limit = 0)
    tss1 = tss[tss[!,"TSS1_r2"].>limit, :]
    tss2 = tss[tss[!,"TSS2_r2"].>limit, :]
    tss2[!,"TSS1_r2"] = tss2[!,"TSS2_r2"]
    tsss = vcat(tss1, tss2)
end


function plot_wt_burstsizes(GENE; limit = 1, r = "r2", maxntoconsider = 50)
    set = sort!(get_genedata(GENE), :Timepoint)
    if GENE == "Hprt" ||GENE == "Fh1"
        set = set[set[!, :Timepoint].!="0i", :]
        set = set[set[!, :Timepoint].!="8i", :]
        set[!, :Timepoint] = [parse(Int, ii)*60 for ii in set[!, :Timepoint]]
    end
    wt_data = set[set[!,:Genotype].=="WT", :]
    wt_data = tss_data(wt_data; limit = limit)
    wt_data = wt_data[wt_data[!,"TSS1_r2"].>limit, :]
    reps = split_by(wt_data, :Rep)
    
    wt_means = vcat([calculate_medians(reps[ii], "TSS1_"*r; col = "Timepoint", rep = ii) for ii in keys(reps)]...)
    
     reps = split_by(wt_data, :Rep)
    
    
    for ii in  keys(reps)
        samp = split_by(reps[ii], :Sample) 
        for aa in   keys(samp)
            if nrow(samp[aa]) > maxntoconsider
             samp[aa] =  samp[aa][shuffle(1:nrow(samp[aa])), :][1:maxntoconsider, :]
            end
        end
        reps[ii] =  vcat([samp[aa] for aa in   keys(samp)]...)
    end
    
    reps = vcat([reps[aa] for aa in   keys(reps)]...)
   
    
    Seaborn.stripplot(data = Pandas.DataFrame(reps), x = "Timepoint", y = "TSS1_"*r, palette = "Set2", hue = "Rep", dodge = false)
    Seaborn.boxplot(data = Pandas.DataFrame(reps), x = "Timepoint", y = "TSS1_"*r, palette = "Blues", showfliers = false)
   
    pretty_axes2()
    
    title(GENE)
    ylabel("Burst Size ("*r*")")
    
    xlabel("Time after LPS (min)")
    legend_removal()
    
    return reps

end

function add_tests_bs(genedata, compsplot, hi; u = 3)
    
t = R"""
tb <- $genedata
tb$Timepoint <- as.factor(tb$Timepoint)
tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Timepoint +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Timepoint
    write.csv(a, "BS_test.csv")
"""
t = CSV.read("BS_test.csv", DataFrames.DataFrame)
    
    pvals =  vcat([t[t[!,:Column1].==ii, :] for ii in compsplot]...)[!,"p adj"]

for ii in 1:length(pvals)
    h = hi[ii]
    plot([ii-1, ii], [h, h], c = "Black")
    annotate(NoLongerProblems.transform_pvalue_in_stars(pvals[ii]), xy =[mean([ii-1, ii]), h+u], ha = "center", va = "center")
    
end
    end


function plot_wt_burstsizes(df; limit = 0, r = "r2", maxntoconsider = 50, GENE = "", max = 10000)
    set = df
    if GENE == "Hprt" ||GENE == "Fh1"
        set = set[set[!, :Timepoint].!="0i", :]
        set = set[set[!, :Timepoint].!="8i", :]
        set[!, :Timepoint] = [parse(Int, ii)*60 for ii in set[!, :Timepoint]]
    end
    wt_data = set
    wt_data = tss_data(wt_data; limit = limit)
    wt_data = wt_data[wt_data[!,"TSS1_"*r].>limit, :]
    reps = split_by(wt_data, :Rep)
    
    wt_means = vcat([calculate_medians(reps[ii], "TSS1_"*r; col = "Timepoint", rep = ii) for ii in keys(reps)]...)
    
     reps = split_by(wt_data, :Rep)
    
    
    for ii in  keys(reps)
        samp = split_by(reps[ii], :Sample) 
        for aa in   keys(samp)
            if nrow(samp[aa]) > maxntoconsider
             samp[aa] =  samp[aa][shuffle(1:nrow(samp[aa])), :][1:maxntoconsider, :]
            end
        end
        reps[ii] =  vcat([samp[aa] for aa in   keys(samp)]...)
    end
    
    reps = vcat([reps[aa] for aa in   keys(reps)]...)
   
    
    Seaborn.stripplot(data = Pandas.DataFrame(reps), x = "Timepoint", y = "TSS1_"*r, hue = "Genotype", dodge = true)
    Seaborn.boxplot(data = Pandas.DataFrame(reps), x = "Timepoint", y = "TSS1_"*r, hue = "Genotype", showfliers = false, dodge = true)
   
    pretty_axes2()
    
    title(GENE)
    ylabel("Burst Size ("*r*")")
    
    xlabel("Time after LPS (min)")
    legend_removal()
    ylim(0, max)
    
    return reps

end

function add_tests_bs(genedata, compsplot, hi; u = 3)
    
t = R"""
tb <- $genedata
tb$Timepoint <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Timepoint +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Timepoint
    write.csv(a, "BS_test.csv")
"""
t = CSV.read("BS_test.csv", DataFrames.DataFrame)
    
    pvals =  vcat([t[t[!,:Column1].==ii, :] for ii in compsplot]...)[!,"p adj"]
    

cords = [
        [-0.25, 0.25],
        [0.75, 1.25],
        [-0.25, 0.75],
        [0.25, 1.25]
    ]
for ii in 1:length(pvals)
    h = hi[ii]
    plot(cords[ii], [h, h], c = "Black")
    annotate(NoLongerProblems.transform_pvalue_in_stars(pvals[ii]), xy =[mean(cords[ii]), h+u], ha = "center", va = "center")
    
end
end


function BS_test(genedata; gene = "Egr2_Ehn")
    
t = R"""
tb <- $genedata
tb$Timepoint <- as.factor(tb$Timepoint)
tb$Sample <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Sample +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Sample
    write.csv(a, paste0($AnalysisName, "/", $gene, "_BS_test.csv"))
"""
t = CSV.read("$AnalysisName"* "/"* "$gene"* "_BS_test.csv", DataFrames.DataFrame)
    
end


function do_mantelhaen(df, comp1, comp2; col = :N_TSS)
    s1 = df[df[!,:Sample] .== comp1, :]
    s2 = df[df[!,:Sample] .== comp2, :]
    tb_test = join_in_all_common_columns(s1, s2)
    tb1 = tb_test; tb1[!,:Burst] = ["Yes" for ii in 1:nrow(tb1)]
    tb1[!,:Count] = tb1[!,col]
    tb1 = tb1[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb2 = tb_test; tb2[!,:Burst] = ["No" for ii in 1:nrow(tb1)]
    tb2[!,:Count] = 2tb2[!,:N_cells] .- tb2[!,col]
    tb2 = tb2[:, [:Sample, :Count, :Rep, :Burst]]
    
    tb = join_in_all_common_columns(tb1, tb2)

    
test = R"""
#library("psych")
#library("vcd")
#library("DescTools")
#library("rcompanion")
library("dplyr")

tb = $tb
tb$Count = as.integer(tb$Count) 

Data <- mutate(tb,
           Sample = factor(Sample, levels=unique(Sample)),
           Burst = factor(Burst, levels=unique(Burst)),
           Rep = factor(Rep, levels=unique(Rep))
           )

# Last variable is the strata (the variable that is not check for assotiation)
Data.xtabs <- xtabs(Count ~ Burst + Sample + Rep, 
                       data=Data)

ftable(Data.xtabs)  

mantelhaen.test(Data.xtabs)
"""

return test[3][1],test[5][1]
end


function BF_test(df, gene, col, comps)
     al = df


odds = []
p = []

for comp in comps
    comp1 = comp[1]
    comp2 = comp[2]
    
    pval = do_mantelhaen(al, comp1, comp2, col = col)[1]; push!(p, pval)
    od = do_mantelhaen(al, comp1, comp2, col = col)[2]; push!(odds, od)
end
    
d = DataFrames.DataFrame()
d[!,:s1] = [comp[1] for comp in comps]
d[!,:s2] = [comp[2] for comp in comps]
d[!,:pval] = p
d[!,:commonoddsratio] = odds
d[!,:padj] = MultipleTesting.adjust([ii for ii in p], MultipleTesting.Bonferroni())

CSV.write("$AnalysisName/"*gene*"_BF_test.csv", d)
d
end

function NoLongerProblems.transform_pvalue_in_stars(p)
    return round(round(p, sigdigits = 2), sigdigits = 3)
end

    
    function add_tests_bs(genedata, compsplot, hi; u = 3)
    
t = R"""
tb <- $genedata
tb$Timepoint <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Timepoint +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Timepoint
    write.csv(a, "BS_test.csv")
"""
t = CSV.read("BS_test.csv", DataFrames.DataFrame)
    
    pvals =  vcat([t[t[!,:Column1].==ii, :] for ii in compsplot]...)[!,"p adj"]
    

cords = [
        [0.75, 1.25],
        [-0.25, 0.75],
        [0.25, 1.25]
    ]
for ii in 1:length(pvals)
    h = hi[ii]
    plot(cords[ii], [h, h], c = "Black")
    annotate(NoLongerProblems.transform_pvalue_in_stars(pvals[ii]), xy =[mean(cords[ii]), h+u], ha = "center", va = "center")
    
end
end

function add_test(df, s1, s2, h, u, cords)
    p = 1
    if sum(t[!,:Column1].== s1*"-"*s2)==1
    p = t[t[!,:Column1].== s1*"-"*s2, "p adj"]
    else
    p = t[t[!,:Column1].== s2*"-"*s1, "p adj"]
    end
    plot(cords, [h, h], c = "Black")
    annotate(NoLongerProblems.transform_pvalue_in_stars(p[1]), xy =[mean(cords), h+u], ha = "center", va = "center")

end
        
        
function plot_burstsizes(df; limit = 0, r = "r2", maxntoconsider = 50, GENE = "", max = 150,  palette = ["red", "darkgray"], order = ["DMSO", "DMSO-LPS","BD1-LPS","BD2-LPS"])
    set = df
    if GENE == "Hprt" ||GENE == "Fh1"
        set = set[set[!, :Timepoint].!="0i", :]
        set = set[set[!, :Timepoint].!="8i", :]
        set[!, :Timepoint] = [parse(Int, ii)*60 for ii in set[!, :Timepoint]]
    end
    wt_data = set
    wt_data = tss_data(wt_data; limit = limit)
    wt_data = wt_data[wt_data[!,"TSS1_"*r].>limit, :]
    reps = split_by(wt_data, :Rep)
    
    wt_means = vcat([calculate_medians(reps[ii], "TSS1_"*r; col = "Timepoint", rep = ii) for ii in keys(reps)]...)
    
     reps = split_by(wt_data, :Rep)
    
    
    for ii in  keys(reps)
        samp = split_by(reps[ii], :Sample) 
        for aa in   keys(samp)
            if nrow(samp[aa]) > maxntoconsider
             samp[aa] =  samp[aa][shuffle(1:nrow(samp[aa])), :][1:maxntoconsider, :]
            end
        end
        reps[ii] =  vcat([samp[aa] for aa in   keys(samp)]...)
    end
    
    reps = vcat([reps[aa] for aa in   keys(reps)]...)
   
    
    #Seaborn.stripplot(data = Pandas.DataFrame(reps), x = "Timepoint", y = "TSS1_"*r, hue = "Rep", palette = "Greys", order = order, jitter = 0.5, zorder = 0)
    Seaborn.boxplot(data = Pandas.DataFrame(reps), x = "Timepoint", y = "TSS1_"*r, hue = "Genotype", showfliers = false, palette = palette,dodge = true, order = order)
   
    pretty_axes2();line075black();squareplot();

    
    title(GENE)
    ylabel("Burst Size ("*r*")")
    
    xlabel("Time after LPS (min)")
    legend_removal()
    ylim(0, max)
    
    return reps


end

function BF_as_percent_max(tb)
    sptb = split_by(tb, :Rep)
    
    dfs = []
    for ii in keys(sptb)
        sptb[ii][!,:BF_percentmax] = sptb[ii][!,:BF] .- sptb[ii][sptb[ii].=="WT_0",:BF_percentmax][1]
        sptb[ii][!,:BF_percentmax] = sptb[ii][!,:BF_percentmax] ./ (maximum(sptb[ii][!,:BF_percentmax])- sptb[ii][sptb[ii].=="WT_0",:BF_percentmax][1])
  
    end
    return   vcat([sptb[ii] for ii in keys(sptb)]...)
end





function BF_as_percent_max(tb, samplemax, samplemin)
    sptb = split_by(tb, :Sample)
    maximum = mean(sptb[samplemax][!,:BF])
    minimum = mean(sptb[samplemin][!,:BF])
    range = maximum - minimum
    dfs = []
    tb[!,:BF_percentmax] = (tb[!,:BF] .-  minimum) ./ range .*100
   return tb
end

        
function plot_BF_enhancer_gene_as_percent_max(tbgene, tbenhancer; genename = "Gene", enhname = "Enhancer", ci = 70)
    genedata = BF_as_percent_max(tbgene, "WT_90", "WT_0")
    enhdata = BF_as_percent_max(tbenhancer, "WT_60", "WT_0")

    genedata[!,:Target] = [genename for ii in 1:nrow(genedata)]
    enhdata[!,:Target] = [enhname for ii in 1:nrow(enhdata)]
    tb  = vcat(genedata, enhdata)
    
    pd = Pandas.DataFrame(tb)
    
    py"""
    
    import seaborn
    
    seaborn.lineplot(x="Timepoint", y="BF_percentmax",
             hue="Target",markers=True, ci = $ci,
             data= $pd)
    

    
    """
    
    ylabel("% of burst frequency increase")
    xlabel("Time after LPS (min)")
    legend_out_of_plot();pretty_axes2();line075black();squareplot();xticks(rotation=0)
    

   print( R"""
    aov.s = aov(BF_percentmax ~ Target*Timepoint  ,data=$tb)
    summary(aov.s)
    """)
    
    #CSV.read(filename, DataFrames.DataFrame)
end
        
        