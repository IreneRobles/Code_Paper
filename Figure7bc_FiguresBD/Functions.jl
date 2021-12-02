
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

function calculatesumary_byrep(exps; r = 2, limit = 1)
expe =  [exps[exps[!,:Rep].==ii, :] for ii in unique(exps[!,:Rep])]
   join_in_all_common_columns([calculatesumary(ex; r = r, limit = limit)  for ex in expe]...)
end


function do_mantelhaen(df, comp1, comp2; col = :N_TSS)
    s1 = df[df[!,:Sample] .== comp1, :]
    s2 = df[df[!,:Sample] .== comp2, :]
    tb_test = join_in_all_common_columns(s1, s2)
    tb1 = tb_test; tb1[!,:Burst] = ["Yes" for ii in 1:nrow(tb1)]
    tb1[!,:Count] =  tb1[!,col]
    tb1 = tb1[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb2 = tb_test; tb2[!,:Burst] = ["No" for ii in 1:nrow(tb1)]
    tb2[!,:Count] = 2tb2[!,:N_cells] .- tb2[!,col]
    tb2 = tb2[:, [:Sample, :Count, :Rep, :Burst]]
    
    tb = join_in_all_common_columns(tb1, tb2)

    try
test = R"""
library("dplyr")
options(warn=-1)
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
    catch
        return 1.,1.
    end

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

CSV.write(gene*"_BF_test.csv", d)
d
end

function NoLongerProblems.transform_pvalue_in_stars(p)
    return round(p, sigdigits = 3)
end
    

function tss_data(tss; limit = 0)
    tss1 = tss[tss[!,"TSS1_r2"].>limit, :]
    tss2 = tss[tss[!,"TSS2_r2"].>limit, :]
    tss2[!,"TSS1_r2"] = tss2[!,"TSS2_r2"]
    tsss = vcat(tss1, tss2)
end


    
function plot_burstsizes_genotype(df, genotype; limit = 1, r = "r2", maxntoconsider = 50, GENE = "", ylimit = 150,  palette = [ "darkgray","red"], order = ["DMSO", "DMSO-LPS","BD1-LPS","BD2-LPS"])
    set = df[df[!,:Genotype].== genotype, :]
    tss_data_ = tss_data(set; limit = limit)
    Seaborn.boxplot(data = Pandas.DataFrame(tss_data_), x = "Timepoint", y = "TSS1_"*r, hue = "Genotype", showfliers = false, palette = palette,dodge = true, order = order)
    title(GENE)
    ylabel("Burst Size ("*r*")")
    xlabel("Treatment")
    legend_removal()
    ylim(0, ylimit)
    pretty_axes2();line075black();squareplot();xticks(rotation=90)
end

function findpadj(t,s1,s2; xy = (0,1), h = 60, u = 5)
p = 1
if sum(t[!,:Column1].== s1*"-"*s2)==1
    p = t[t[!,:Column1].== s1*"-"*s2, "p adj"][1]
else
    p = t[t[!,:Column1].== s2*"-"*s1, "p adj"][1]
end
        xy2 = [Statistics.mean(xy), h+u]
        plt.plot(xy, [h, h], c = "black", lw = 0.75)

        plt.annotate(NoLongerProblems.transform_pvalue_in_stars(p), xy = xy2, ha= "center")
        return p
        
end
    
function plot_BF_Genotype(BFs, TEST; ylimit = 0.7, 
        hs = [0.25, 0.55, 0.4, 0.50],
        ylabel_ = "BF", Genotype = "WT",
        colorstripplot = "gray",
        colorboxplot = "darkgray",
    )
    BFssub = BFs[BFs[!,:Genotype].==Genotype, :]
    
Seaborn.stripplot(data = Pandas.DataFrame(BFssub), y = "BF", x = "Timepoint", hue = "Genotype", dodge = true, palette = [colorstripplot, "red"],
    order = ["DMSO", "DMSO-LPS","BD1-LPS","BD2-LPS"])
Seaborn.boxplot(data = Pandas.DataFrame(BFssub), y = "BF", x = "Timepoint", hue = "Genotype", dodge = true, palette = [colorboxplot, "red"],
    order = ["DMSO", "DMSO-LPS","BD1-LPS","BD2-LPS"],showfliers = false)
    xlabel("Treatment")

    
    xys =[ [0, 1], [1,2], [1, 3], [2,3]]
    
    o = if Genotype == "WT" 5 else 9 end
    u = ylimit/30
    for ii in 0:(length(hs)-1)
        h = hs[ii+1]
      
        xy =xys[ii+1] ; xy2 = [Statistics.mean(xy), h+u]
        plt.plot(xy, [h, h], c = "black")
        pval = TEST[ii+o, :padj]
        plt.annotate(NoLongerProblems.transform_pvalue_in_stars(pval), xy = xy2, ha= "center")
    end
    ylabel(ylabel_)
    ylim(0, ylimit);legend_removal();pretty_axes2();line075black();squareplot();xticks(rotation=90)

TEST
end
    