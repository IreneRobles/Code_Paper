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
using Statistics
using PyPlot
using PrettyPlotting
using DataFrames
include("../TSS_quantification/StandardCode.jl")
using RCall
using Distributions


function calculate_bf(t, tname; limit = 2)
    new_df = DataFrames.DataFrame()
    
    samples = unique(t[!,:Sample])
    new_df[!,:Sample] = samples
    new_df[!,:Genotype] = [split(ii, "_")[1] for ii in samples]
    new_df[!,:Timepoint] = [split(ii, "_")[2] for ii in samples]
    
    n_cells = []
    cellsactive1 = []
    cellsactive2 = []
    n_tss =[]
    bs_mean =[]
    bs_median =[]
     bs_std =[]
    
    for ii in samples
        push!(n_cells, sum(t[!,:Sample] .== ii))
        sp_sam = t[t[!,:Sample] .== ii, :]
        locus1_bool = sp_sam[!,Symbol(string("locus1_", tname, "_size"))] .> limit
        locus1 = sp_sam[locus1_bool,Symbol(string("locus1_", tname, "_size"))]
        locus2_bool = sp_sam[!,Symbol(string("locus2_", tname, "_size"))] .> limit
        locus2 = sp_sam[locus2_bool,Symbol(string("locus2_", tname, "_size"))]
        push!(cellsactive1, Statistics.mean(locus1_bool.|locus2_bool))
        push!(cellsactive2, Statistics.mean(locus1_bool.*locus2_bool))
        push!(n_tss, length(locus1)+length(locus2))
        sizes = append!(locus1, locus2)
        push!(bs_mean, Statistics.mean(sizes))
        mediana = if isempty(sizes) 0 else median(sizes) end
        push!(bs_median, mediana)
        push!(bs_std, Statistics.std(sizes))
        
    end
    
    new_df[!,:N_cells] = n_cells
    new_df[!,Symbol(string("CellsActive1_", tname))] = cellsactive1
    new_df[!,Symbol(string("CellsActive2_", tname))] = cellsactive2
    new_df[!,Symbol(string(tname, "_N"))] = n_tss
    
    new_df[!,Symbol(string("BF_", tname))] = n_tss./2n_cells
    new_df[!,Symbol(string("BS_mean_", tname))] = bs_mean
    new_df[!,Symbol(string("BS_median_", tname))] = bs_median
    new_df[!,Symbol(string("BS_std_", tname))] = bs_std
    
    new_df
    
end



function calculate_bf_by_rep(df, tname; limit = 2)
    
    df[!,:Sample] = [string(df[ii, :Genotype], "_",df[ii, :Timepoint]) for ii in 1:nrow(df)]
    subdfs = NoLongerProblems.split_by(df, :Rep)
    reps = collect(keys(subdfs))
    
    subdfs = [calculate_bf(subdfs[ii], tname; limit = limit) for ii in reps]
    
    
    for ii in 1:length(reps)
          subdfs[ii][:Rep]= reps[ii]
    end
    
    join_in_all_common_columns(subdfs)
end



function calculate_randoms(coocurrences, probe1, probe2)
    
    coocurrences[!,Symbol("ByCell_"*probe1*"_"*probe2*"_random")] = coocurrences[!,"CellsAct_"*probe1].* coocurrences[!,"CellsAct_"*probe2]
    coocurrences[!,Symbol("ByLocus_"*probe1*"_"*probe2*"_random")] = coocurrences[!,"BF_"*probe1].* coocurrences[!,"BF_"*probe2]

    coocurrences
    
end






function get_locus_data(exp, tss)
    cols1 = NoLongerProblems.columns_containing(exp, "locus1_"*lowercase(tss))
    cols11 = NoLongerProblems.columns_containing(exp, "locus1_"*uppercase(tss))
    cols1 = intersect(cols1, cols11)
    cols1 = push!([Symbol(ii) for ii in cols1], :Cell, :Sample, :Genotype, :Timepoint, :Image, :Rep, :Sample_Rep)
    locus1 = exp[:, cols1]
    cols2 = NoLongerProblems.columns_containing(exp, "locus2_"*lowercase(tss))
    cols21 = NoLongerProblems.columns_containing(exp, "locus2_"*uppercase(tss))    
    cols2 = intersect(cols2, cols21)
    cols2 = push!([Symbol(ii) for ii in cols2], :Cell, :Sample, :Genotype, :Timepoint, :Image, :Rep, :Sample_Rep)
    locus2 = exp[:, cols2]
    for col in cols2
        rename!(locus2, col => Symbol(string(replace(string(col), "locus2" => "locus1"))))
    end
    return join_in_all_common_columns(locus1, locus2)
end


function do_mantelhaen(df, comp1, comp2; col = :N_TSS2)
    s1 = df[df[!,:Sample] .== comp1, :]
    s2 = df[df[!,:Sample] .== comp2, :]
    tb_test = join_in_all_common_columns(s1, s2)
    tb1 = tb_test; tb1[:Burst] = ["Yes" for ii in 1:nrow(tb1)]
    tb1[!,:Count] = tb1[!,col]
    tb1 = tb1[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb2 = tb_test; tb2[:Burst] = ["No" for ii in 1:nrow(tb1)]
    tb2[!,:Count] = tb2[!,:N_cells] * 2 - tb2[!,col]
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


function BF_test(df, gene, col )
    al = df

comps = [
    ["WT_0","WT_120"],
    ["Rad21KO_0","Rad21KO_120"],
    ["WT_0","Rad21KO_0"],
    ["WT_120", "Rad21KO_120"],

]


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

CSV.write("BF_fig/"*gene*"_BF_test.csv")
d
end


function BS_test_TSS2(locus; gene = "gene")

CSV.write("temp.csv", locus)


t = R"""


tb <- read.csv("temp.csv")

tb$Genotype <- as.factor(tb$Genotype)
tb$Timepoint <- as.factor(tb$Timepoint)
tb$Rep <- as.factor(tb$Rep);

s <- tb
aov.s = aov(locus1_TSS2_size ~ Sample + Rep  ,data=s)  #do the analysis of variance
test <- TukeyHSD(aov.s)
test$Sample
write.csv(test$Genotype,"TukeyHSD_Fh1.csv")
a = test$Sample
b = $gene
names(a)[names(a) == "p adj"] <- "padj"
write.csv(a, paste0("BS_test/", b,"_BS_test.csv"))

"""
t = DataFrame!(CSV.File("BS_test/"*gene*"_BS_test.csv"))

rename!(t, Symbol("p adj")=> :padj)
end

function BS_test_TSS4(locus; gene = "gene")

CSV.write("temp.csv", locus)


t = R"""


tb <- read.csv("temp.csv")

tb$Genotype <- as.factor(tb$Genotype)
tb$Timepoint <- as.factor(tb$Timepoint)
tb$Rep <- as.factor(tb$Rep);

s <- tb
aov.s = aov(locus1_TSS4_size ~ Sample + Rep  ,data=s)  #do the analysis of variance
test <- TukeyHSD(aov.s)
test$Sample
write.csv(test$Genotype,"TukeyHSD_Fh1.csv")
a = test$Sample
b = $gene
names(a)[names(a) == "p adj"] <- "padj"
write.csv(a, paste0("BS_test/", b,"_BS_test.csv"))

"""
t = DataFrame!(CSV.File("BS_test/"*gene*"_BS_test.csv"))

rename!(t, Symbol("p adj")=> :padj)
end
function BS_test_TSS3(locus; gene = "gene")

CSV.write("temp.csv", locus)


t = R"""


tb <- read.csv("temp.csv")

tb$Genotype <- as.factor(tb$Genotype)
tb$Timepoint <- as.factor(tb$Timepoint)
tb$Rep <- as.factor(tb$Rep);

s <- tb
aov.s = aov(locus1_TSS3_size ~ Sample + Rep  ,data=s)  #do the analysis of variance
test <- TukeyHSD(aov.s)
test$Sample
write.csv(test$Genotype,"TukeyHSD_Fh1.csv")
a = test$Sample
b = $gene
names(a)[names(a) == "p adj"] <- "padj"
write.csv(a, paste0("BS_test/", b,"_BS_test.csv"))

"""
t = DataFrame!(CSV.File("BS_test/"*gene*"_BS_test.csv"))

rename!(t, Symbol("p adj")=> :padj)
end

function calculate_allbf(df;limit = 1)
    genes = unique([split(split(ii, "locus")[end][3:end], "_size")[1] for ii in NoLongerProblems.columns_containing(df, "_size")])
    datatoreturn = [sort!(calculate_bf_by_rep(df, gene; limit = limit), [:Sample, :Rep]) for gene in genes]
    return datatoreturn 
end

function merge_BF(listbfs)
    
    for bf in listbfs
        bf[!,:Sample_Rep] = bf[!,:Sample] .* "_" .* string.(bf[!,:Rep] )
    end
    
    bfs = listbfs[1]
    
    for d in listbfs[2:end]
        cols = [:Sample_Rep, Symbol.(NoLongerProblems.columns_containing(d, "TSS"))...]
        d=d[!,cols]
        bfs = join(bfs, d, on = :Sample_Rep)
    
    end
    
    bfs
    
end

function calculate_coocurrences(DF, tss1, tss2; limit1 = 1, limit2 = 0, col = :Sample_Rep)
    DF[!,:Sample_Rep] = DF[!,:Sample].*"_".*string.(DF[!,:Rep])
    sams = unique(DF[!,col])
     
    bycell =[] 
    bylocus = []
    tss1s = []
    tss2s = []
    cell1 = []
    cell2 = []
     ncells = []
    
    
    for ii in sams
        subDF = DF[DF[!,col] .== ii, :]
        
        locus1_bool_tss1 = subDF[!,Symbol("locus1_", tss1, "_size")].> limit1
        locus1_bool_tss2 = subDF[!,Symbol("locus1_", tss2, "_size")].> limit2
        
        locus2_bool_tss1 = subDF[!,Symbol("locus2_", tss1, "_size")].> limit1
        locus2_bool_tss2 = subDF[!,Symbol("locus2_", tss2, "_size")].> limit2
        
        locus1_co_bool = locus1_bool_tss1 .* locus1_bool_tss2
        locus2_co_bool = locus2_bool_tss1 .* locus2_bool_tss2
        
        locus_co_avg = Statistics.mean(append!(locus1_co_bool,locus2_co_bool))
        push!(bylocus, locus_co_avg)
        

        
        tss1_bool = locus1_bool_tss1 .| locus2_bool_tss1
        tss2_bool = locus1_bool_tss2 .| locus2_bool_tss2
        
        
        bycell_co = Statistics.mean(tss1_bool.*tss2_bool)
        
        push!(bycell, bycell_co)
        push!(tss1s, Statistics.mean(append!(locus1_bool_tss1, locus2_bool_tss1)))
        push!(tss2s, Statistics.mean(append!(locus1_bool_tss2, locus2_bool_tss2)))
        push!(cell1, Statistics.mean(tss1_bool))
        push!(cell2, Statistics.mean(tss2_bool))
        push!(ncells, length(tss2_bool))
       
        
    end
    
    newdf = DataFrames.DataFrame(Dict(
         col => sams,
         Symbol("ByLocus_"*tss1*"_"*tss2)   => bylocus,
            Symbol("BF_"*tss1)=> tss1s,
            Symbol("BF_"*tss2)=> tss2s,
            Symbol("CellsAct_"*tss1)=> cell1,
            Symbol("CellsAct_"*tss2)=> cell2,
        Symbol("ByCell_"*tss1*"_"*tss2)   => bycell,
            Symbol("N_cells")=>ncells
   

        )
    )
    
    
end



function calculate_coocurrences(DF, tss1, tss2, tss3; limit = 2, col = :Sample_Rep)
    DF[!,:Sample_Rep] = DF[!,:Sample].*"_".*string.(DF[!,:Rep])
    sams = unique(DF[!,col])
    
    bycell =[] 
    bylocus = []
    
    for ii in sams
        subDF = DF[DF[!,:Sample_Rep] .== ii, :]
        
        locus1_bool_tss1 = subDF[!,Symbol("locus1_", tss1, "_size")].> limit
        locus1_bool_tss2 = subDF[!,Symbol("locus1_", tss2, "_size")].> limit
        locus1_bool_tss3 = subDF[!,Symbol("locus1_", tss3, "_size")].> limit
        
        locus2_bool_tss1 = subDF[!,Symbol("locus2_", tss1, "_size")].> limit
        locus2_bool_tss2 = subDF[!,Symbol("locus2_", tss2, "_size")].> limit
        locus2_bool_tss3 = subDF[!,Symbol("locus2_", tss3, "_size")].> limit
        
        locus1_co_bool = locus1_bool_tss1 .* locus1_bool_tss2 .* locus1_bool_tss3
        locus2_co_bool = locus2_bool_tss1 .* locus2_bool_tss2 .* locus2_bool_tss3
        
        locus_co_avg = Statistics.mean(append!(locus1_co_bool,locus2_co_bool))
        push!(bylocus, locus_co_avg)
        
        tss1_bool = locus1_bool_tss1 .| locus2_bool_tss1
        tss2_bool = locus1_bool_tss2 .| locus2_bool_tss2
        tss3_bool = locus1_bool_tss3 .| locus2_bool_tss3
        
        bycell_co = Statistics.mean(tss1_bool.*tss2_bool.*tss3_bool)
        
        push!(bycell, bycell_co)

        
    end
    
    newdf = DataFrames.DataFrame(Dict(
         col => sams,
         Symbol("ByLocus_"*tss1*"_"*tss2*"_"*tss3)   => bylocus,
        Symbol("ByCell_"*tss1*"_"*tss2*"_"*tss3)   => bycell

        )
    )
    
    
end




function Rcortest(xi, yi)
    test = R"cor.test($xi, $yi)"
end


function bernoulli_prob(ncells, bf1, bf2; iterations = 100000)
    #bf1 and bf2 correspond to probabilities between 0 and 1 (fraction of actie cells or burst fraction)
    
    coocs = Array{Float64,1}()
    
    for ii in 1:iterations
    
        bernu1 = Bernoulli(bf1)

        # Simulate sampling from gene1
        gene1 = rand(bernu1, ncells)
        
        # Simulate sampling from gene2
        
        bernu2 = Bernoulli(bf2)
        gene2 = rand(bernu2, ncells)
        
        # Calculate coocurrence
        
        cooc = mean(gene1 .* gene2)
        
        append!(coocs, cooc)
        
    end
    
    coocs
    
end

function bernoulli_prob(ncells, bf1, bf2, bf3; iterations = 100000)
    #bf1 and bf2 correspond to probabilities between 0 and 1 (fraction of actie cells or burst fraction)
    
    coocs = Array{Float64,1}()
    
    for ii in 1:iterations
    
        bernu1 = Bernoulli(bf1)
        # Simulate sampling from gene1
        gene1 = rand(bernu1, ncells)
        
        # Simulate sampling from gene2
        
        bernu2 = Bernoulli(bf2)
        gene2 = rand(bernu2, ncells)
        
        
         bernu3 = Bernoulli(bf3)
        # Simulate sampling from gene1
        gene3 = rand(bernu3, ncells)
        
        # Calculate coocurrence
        
        cooc = mean(gene1 .* gene2.* gene3)
        
        append!(coocs, cooc)
        
    end
    
    coocs
    
end


function calculate_pvalue_cooc_mut_ex(value, randomsampled; iterations = 100000, plotthings = false, c = "red", label = "Exp", dist = true)
    
    if plotthings
        if dist 
        plt.hist(randomsampled, 20, normed = true, color ="lightgray")
        end

        plt.axvline(value,  color = c, label = label)

        plt.legend()
        plt.show()
    end
    
    
    pval_cooc = Statistics.mean(randomsampled .> value)
    pval_mutex = 1 - pval_cooc
    
    
    return pval_cooc, pval_mutex
 
end

function calculate_empirical_pvalues(df, tss1, tss2; iterations = 10000)
    bycellcol = Symbol("ByCell_"*tss1*"_"*tss2)
    bylocuscol = Symbol("ByLocus_"*tss1*"_"*tss2)
    
    pval_cooc_cell = []
    pval_mutexc_cell = []
    
    pval_cooc_loc = []
    pval_mutexc_loc = []
    
    p = ProgressMeter.Progress(nrow(df)*2)
    
    for ii in 1:nrow(df)
        
        
        cell1 = df[ii,Symbol("CellsAct_"*tss1)]
        cell2 = df[ii,Symbol("CellsAct_"*tss2)]
        loc1 = df[ii,Symbol("BF_"*tss1)]
        loc2 = df[ii,Symbol("BF_"*tss2)]
        bycellexp = df[ii,bycellcol]
        bylocexp = df[ii,bylocuscol]
        ncells = df[ii,Symbol("N_cells")]
        
        testcells1 = bernoulli_prob(ncells, cell1, cell2, iterations = iterations)
        testcells1 = calculate_pvalue_cooc_mut_ex(bycellexp, testcells1; iterations = iterations, plotthings = false)
        push!(pval_cooc_cell, testcells1[1]); push!(pval_mutexc_cell, testcells1[2]);
        next!(p)
        
                
        testcells2 = bernoulli_prob(ncells*2, loc1, loc2, iterations = iterations)
        testcells2 = calculate_pvalue_cooc_mut_ex(bylocexp, testcells2; iterations = iterations, plotthings = false)
        push!(pval_cooc_loc, testcells2[1]); push!(pval_mutexc_loc, testcells2[2]);
        next!(p)
        
    end
    
    df[!,Symbol("ByCell_"*tss1*"_"*tss2*"_pval_cooc")] = pval_cooc_cell
    df[!,Symbol("ByCell_"*tss1*"_"*tss2*"_pval_mutexc")] = pval_mutexc_cell
    
    df[!,Symbol("ByLocus_"*tss1*"_"*tss2*"_pval_cooc")] = pval_cooc_loc
    df[!,Symbol("ByLocus_"*tss1*"_"*tss2*"_pval_mutexc")] = pval_mutexc_loc
    
    return df
end

function calculate_empirical_pvalues(df, tss1, tss2, tss3; iterations = 10000)
    bycellcol = Symbol("ByCell_"*tss1*"_"*tss2*"_"*tss3)
    bylocuscol = Symbol("ByLocus_"*tss1*"_"*tss2*"_"*tss3)
    
    pval_cooc_cell = []
    pval_mutexc_cell = []
    
    pval_cooc_loc = []
    pval_mutexc_loc = []
    
    p = ProgressMeter.Progress(nrow(df)*2)
    
    for ii in 1:nrow(df)
        
        
        cell1 = df[ii,Symbol("CellsActive1_"*tss1)]
        cell2 = df[ii,Symbol("CellsActive1_"*tss2)]
        cell3 = df[ii,Symbol("CellsActive1_"*tss3)]
        loc1 = df[ii,Symbol("BF_"*tss1)]
        loc2 = df[ii,Symbol("BF_"*tss2)]
        loc3 = df[ii,Symbol("BF_"*tss3)]
        bycellexp = df[ii,bycellcol]
        bylocexp = df[ii,bylocuscol]
        ncells = df[ii,Symbol("N_cells")]
        
        testcells1 = bernoulli_prob(ncells, cell1, cell2, cell3, iterations = iterations)
        testcells1 = calculate_pvalue_cooc_mut_ex(bycellexp, testcells1; iterations = iterations, plotthings = false)
        push!(pval_cooc_cell, testcells1[1]); push!(pval_mutexc_cell, testcells1[2]);
        next!(p)
        
                
        testcells2 = bernoulli_prob(ncells*2, loc1, loc2, loc3, iterations = iterations)
        testcells2 = calculate_pvalue_cooc_mut_ex(bylocexp, testcells2; iterations = iterations, plotthings = false)
        push!(pval_cooc_loc, testcells2[1]); push!(pval_mutexc_loc, testcells2[2]);
        next!(p)
        
    end
    
    df[!,Symbol("ByCell_"*tss1*"_"*tss2*"_"*tss3*"_pval_cooc")] = pval_cooc_cell
    df[!,Symbol("ByCell_"*tss1*"_"*tss2*"_"*tss3*"_pval_mutexc")] = pval_mutexc_cell
    
    df[!,Symbol("ByLocus_"*tss1*"_"*tss2*"_"*tss3*"_pval_cooc")] = pval_cooc_loc
    df[!,Symbol("ByLocus_"*tss1*"_"*tss2*"_"*tss3*"_pval_mutexc")] = pval_mutexc_loc
    
    return df
end

