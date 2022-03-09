REPS = REPS = [tb[tb[!,:Rep].== ii, :] for ii in [1, 2, 3, 4]]
reps = vcat(REPS...)
reps_expressingexon = sort(reps[reps[!,:N_exon].> lim_exp, :], "Genotype", rev = true)
name = "TukeyHSD_"*gene*"_nexpressing_.csv"
t = R"""
tb <- $reps_expressingexon
tb$Timepoint <- as.factor(tb$Timepoint)
tb$Genotype <- as.factor(tb$Genotype)
tb$Sample <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb
s <- tb
aov.s = aov(N_exon ~ Sample + Rep  ,data=s)  #do the analysis of variance
test <- TukeyHSD(aov.s)
test$Sample
write.csv(test$Sample,$name)
t = test$Sample
"""

t = CSV.read(name, DataFrames.DataFrame)

summary1 = mean_burst_size_and_burst_fraction(REPS..., limit = lim_exp)
summary1 = summary1
summary1[!,:F_Expressing_exon] = summary1[!,:N_Expressing_exon] ./summary1[!,:n_cells]
CSV.write(gene*"_summary_WTRad21KO.csv", summary1)



function do_mantelhaen(df, comp1, comp2)
    s1 = df[df[!,:Sample] .== comp1, :]
    s2 = df[df[!,:Sample] .== comp2, :]
    tb_test = join_in_all_common_columns(s1, s2)
    tb1 = tb_test; tb1[!,:Burst] = ["Yes" for ii in 1:nrow(tb_test)]
    tb1[!,:Count] = tb1[!,:N_Expressing_exon]
    tb1 = tb1[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb2 = tb_test; tb2[!,:Burst] = ["No" for ii in 1:nrow(tb_test)]
    tb2[!,:Count] = tb2[!,:n_cells] .- tb2[!,:N_Expressing_exon]
    tb2 = tb2[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb = join_in_all_common_columns(tb1, tb2)

    
test = R"""

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


   