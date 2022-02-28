ENV["Code"] = "../../Code"
for folder in readdir(ENV["Code"]); push!(LOAD_PATH, normpath(ENV["Code"], folder));end

using CSV
using NoLongerProblems_FileHandling, NoLongerProblems_Pandas, NoLongerProblems
using DataFrames, HypothesisTests, MultipleTesting, Random,Distributions
import Pandas
using PyCall, RCall,PyPlot
using PrettyPlotting,Statistics,Seaborn

include("../TSS_quantification/StandardCode.jl")


suff = "1um_noz"
cooc = []

function get_coocs()

coocs = CSV.read(normpath(genefolder, probe1*"__"*probe2*"__coocurrences"*suff*".csv"), DataFrames.DataFrame)
coocs[!,:BF_Enhancer] = coocs[!,"BF_"*probe2]
coocs[!,:BF_Gene] = coocs[!,"BF_"*probe1]
coocs[!,:Gene] = [probe1 for ii in 1:nrow(coocs)]
if probe1 == "Prdm1_intron"
        coocs[!,:Timepoint] = [split(ii, "_")[2] for ii in coocs[!,:Sample_Rep]]
        coocs = coocs[coocs[!,:Timepoint].!= "30", :]
        coocs = coocs[coocs[!,:Timepoint].!= "90", :]
        
    end

push!(cooc, coocs)
end

genefolder = "../CompleteSets/linkedlocus"
probe1= "Il12b_intron"
probe2 = "HSS1"

get_coocs()

probe1 = "Ifnb1forL2"
probe2 = "L2"

get_coocs()


probe1 = "Egr2_intron"
probe2 = "Egr2_enh"

get_coocs()


probe1 = "Prdm1_intron"
probe2 = "Prdm1_enh"

get_coocs()


probe2 = "Enh"
probe1 = "Peli1_intron"


get_coocs()

tb = join_in_all_common_columns(cooc...)
tb[!,:Genotype] = [split(ii, "_")[1] for ii in tb[!,:Sample_Rep]]
tb[!,:Timepoint] = [split(ii, "_")[2] for ii in tb[!,:Sample_Rep]]

tb[!,"log2_BF_Enhancer"] = log2.(tb[!,"BF_Enhancer"])
tb[!,"log2_BF_Gene"] = log2.(tb[!,"BF_Gene"])



function plotgeneenhancerlog2bf(gene; Genotype = "All")
a = tb

if Genotype != "All"
a = a[a[!,:Genotype].==Genotype, :]
end

figure(figsize = [3,3])
a[!,"log2_BF_Enhancer"] = log2.(a[!,"BF_Enhancer"])
a[!,"log2_BF_Gene"] = log2.(a[!,"BF_Gene"]) 
    

a = a[.! isinf.(a[!,"log2_BF_Gene"]), :]
a = a[.! isinf.(a[!,"log2_BF_Enhancer"]), :]    

x = a[!,"log2_BF_Enhancer"]
y = a[!,"log2_BF_Gene"]
    
CSV.write("../SourceData/Fig3b.csv", a)

pd = Pandas.DataFrame(a)
pd["LPS"] = pd["Timepoint"] .> 0

    
    
if Genotype != "All"
py"""
import seaborn as Seaborn
Seaborn.scatterplot(data = $pd, x = "log2_BF_Enhancer",  y = "log2_BF_Gene", hue = "Gene", s = 100, style = "LPS")
Seaborn.regplot(data = $pd, x = "log2_BF_Enhancer",  y = "log2_BF_Gene", scatter = 0, color = "grey", ci =  95 )
"""
else
        py"""
import seaborn as Seaborn
Seaborn.scatterplot(data = $pd, x = "log2_BF_Enhancer",  y = "log2_BF_Gene", hue = "Genotype",palette = ["black", "red"], s = 100, style = "Gene")
b = $pd
wt = b[b["Genotype"]=="WT"]
    
Seaborn.regplot(data = wt, x = "log2_BF_Enhancer",  y = "log2_BF_Gene", scatter = 0, color = "black", ci =  95 )

rad = b[b["Genotype"]!="WT"]
    
Seaborn.regplot(data = rad, x = "log2_BF_Enhancer",  y = "log2_BF_Gene", scatter = 0, color = "red", ci =  95 )

"""
end


TEST = R"""cor.test($x,$y)"""
corr= round(TEST[4][1], digits = 3)
pval= round(TEST[3][1], sigdigits = 3)
annotate("r = $corr \np = $pval", xy = (-9, -2))
title("$gene")
ylim(-12 , 0.5); xlim(-10, 0.5)




ylabel("log2(Gene \n Burst Frequency)")
xlabel("log2(Enhancer \n Burst Frequency)")
    
legend_out_of_plot()
pretty_axes2()
squareplot()
savefigwithtext(gene*"_BFenhancerBFgene.svg")

return TEST
end

