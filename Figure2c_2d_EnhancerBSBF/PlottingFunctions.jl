function plot_bf(gene; scale = "min")

pd = Pandas.DataFrame(BF[gene])
pd = pd[pd["Genotype"] == "WT"]
Seaborn.boxplot(data = pd, y = "BF_TSS", x = "Timepoint", palette = "Blues", showfliers = false, zorder = 1)

pretty_axes2()
xlabel("Time after LPS ("*scale*")")
ylabel("Burst Freq.")
    title(gene)
    legend_removal()
line075black()
    
end
function do_mantelhaen(df, comp1, comp2)
    s1 = df[df[!,:Sample] .== comp1, :]
    s2 = df[df[!,:Sample] .== comp2, :]
    tb_test = join_in_all_common_columns(s1, s2)
    tb1 = tb_test; tb1[!,:Burst] = ["Yes" for ii in 1:nrow(tb_test)]
    tb1[!,:Count] = tb1[!,:TSS_N]
    tb1 = tb1[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb2 = tb_test; tb2[!,:Burst] = ["No" for ii in 1:nrow(tb_test)]
    tb2[!,:Count] = tb2[!,:N_cells] * 2 .- tb2[!,:TSS_N]
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



function add_tests3(tests, hi; u = 0.005)
    
pvals = adjust([t[1] for t in tests], Bonferroni())

for ii in 1:length(pvals)
    h = hi[ii]
    plot([ii-1, ii], [h, h], c = "Black", lw = 0.75)
        
    p = round(round(pvals[ii], sigdigits = 2), sigdigits = 2)
    annotate(("P = $p"), xy =[mean([ii-1, ii]), h+u], ha = "center", va = "center")
    
end
        ax = gca()
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable="box")
end

function add_tests4(tests, hi; u = 0.005)
    
pvals = [t[1] for t in tests]

for ii in 1:length(pvals)
    h = hi[ii]
    plot([ii-1, ii], [h, h], c = "Black", lw = 0.75)
        
     p = round(round(pvals[ii], sigdigits = 2), sigdigits = 2)
    annotate(("P = $p"), xy =[mean([ii-1, ii]), h+u], ha = "center", va = "center")
    
end
        ax = gca()
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable="box")
end


function add_tests_bs(genedata, comps;  u = 3)
    
genedata[!,:Sample] =  genedata[!,:Genotype]  .* "_" .* string.(genedata[!,:Timepoint]) 
sort!(genedata, :Sample)
    
t = R"""
tb <- $genedata
tb$Timepoint <- as.factor(tb$Sample)


tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Sample +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Sample
    write.csv(a, paste0($gene,"BS_test.csv"))
    a
"""
 a = CSV.read(string(gene,"BS_test.csv"), DataFrame)
    
[a[a[!,:Column1].== ii, "p adj"]   for ii in comps]
end


function plot_bs(gene, limit, time1, time2; scale = "min", maxy = 50)
df = tss_data(get_completeset(gene); limit = limit)
df[!,:Timepoint] = [if ii .== "135" "120" else ii end for ii in df[!,:Timepoint]]
bool1 = df[!,:Timepoint] .== time1 
bool2 = df[!,:Timepoint] .== time2
bool = bool1 .| bool2
pd = Pandas.DataFrame(df[bool, :])
pd = pd[pd["Genotype"] == "WT"]

Seaborn.boxplot(data = pd, y = "TSS1_r2",  x = "Timepoint", showfliers = false, palette = "Blues",)
   # Seaborn.stripplot(data = pd, y = "TSS1_r2", hue = "Genotype", x = "Timepoint", dodge = true, zorder = 0, palette = "Blues", jitter = 0.4, s = 1)

pretty_axes2()
xlabel("Time after LPS ("*scale*")")
ylabel("Burst Size")
    title(gene)
    legend_removal()
    ylim(0, maxy)
    ax = gca()
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable="box")
    line075black()

    return df
    
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


function get_genedata(genesymbol)
    genedata = CSV.read("GeneData/"*genesymbol*".csv", DataFrame)
       try
        genedata[!,:TSS1_r2] = [if ii == "NA" 0 else parse(Float64, ii) end for ii in genedata[!,:TSS1_r2]]
    catch
    end
    
       try
        genedata[!,:TSS2_r2] = [if ii == "NA" 0 else parse(Float64, ii) end for ii in genedata[!,:TSS2_r2]]
    catch
    end
    genedata
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





function plot_wt_burstsizes(GENE; limit = 0, r = "r2", maxntoconsider = 50)
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
aov.s = aov(TSS1_r2 ~ Timepoint + Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Timepoint
    write.csv(a, "BS_test.csv")
"""
t = CSV.read("BS_test.csv", DataFrame)
    
    pvals =  vcat([t[t[!,:Column1].==ii, :] for ii in compsplot]...)[!,"p adj"]

for ii in 1:length(pvals)
    h = hi[ii]
    plot([ii-1, ii], [h, h], c = "Black")
    annotate(NoLongerProblems.transform_pvalue_in_stars(pvals[ii]), xy =[mean([ii-1, ii]), h+u], ha = "center", va = "center")
    
end
    end


function tss_data(tss; limit = 0)
    tss1 = tss[tss[!,"TSS1_r2"].>limit, :]
    tss2 = tss[tss[!,"TSS2_r2"].>limit, :]
    tss2[!,"TSS1_r2"] = tss2[!,"TSS2_r2"]
    tsss = vcat(tss1, tss2)
end