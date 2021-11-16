using Seaborn, NoLongerProblems_Pandas, Pandas, NoLongerProblems, PrettyPlotting, DataFrames, PyPlot
using NoLongerProblems, Statistics
using CSV, DataFrames, MultipleTesting
using FQfiles, RCall, PyCall
    
    
tb = CSV.read("../CompleteSets/GeneData/Il12b_intron_original.csv", DataFrames.DataFrame);
reps = tb
reps = reps[reps[!,:Genotype].=="WT", :]

R"""
tb <- $reps
tb$Timepoint <- as.factor(tb$Timepoint)
tb$Genotype <- as.factor(tb$Genotype)
tb$Sample <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb
s <- tb
aov.s = aov(N_exon ~ Sample + Rep  ,data=s)  #do the analysis of variance
test <- TukeyHSD(aov.s)
test$Sample
write.csv(test$Sample,"TukeyHSD_Il12b_mature_WT.csv")
test$Sample
"""

tb = CSV.read("../CompleteSets/GeneData/Il12b_intron_original.csv", DataFrames.DataFrame);
reps = tb
reps = reps[reps[!,:Genotype].=="Rad21KO", :]

R"""
tb <- $reps
tb$Timepoint <- as.factor(tb$Timepoint)
tb$Genotype <- as.factor(tb$Genotype)
tb$Sample <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb
s <- tb
aov.s = aov(N_exon ~ Sample + Rep  ,data=s)  #do the analysis of variance
test <- TukeyHSD(aov.s)
test$Sample
write.csv(test$Sample,"TukeyHSD_Il12b_mature_Rad21KO.csv")
test$Sample
"""