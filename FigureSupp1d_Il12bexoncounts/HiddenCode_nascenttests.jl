tb = CSV.read("../CompleteSets/GeneData/Il12b_intron_original.csv", DataFrames.DataFrame);

REPS = [tb[tb[!,:Rep].== ii, :] for ii in [1, 2, 3, 4]]

locusd = locusdata2(REPS...)
locusd = locusd[locusd[!,:Genotype].=="WT", :]
limTSS = 1
locusd = locusd[locusd[!,:BurstSize].> limTSS, :]
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
write.csv(test$Sample,"TukeyHSD_Il12b_BurstSize_WT.csv")
test$Sample
"""

locusd = locusdata2(REPS...)
locusd = locusd[locusd[!,:Genotype].=="Rad21KO", :]
limTSS = 1
locusd = locusd[locusd[!,:BurstSize].> limTSS, :]

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
write.csv(test$Sample,"TukeyHSD_Il12b_BurstSize_Rad21KO.csv")
test$Sample
"""

summary1 = mean_burst_size_and_burst_fraction(REPS...)
summary1 = summary1[summary1[!,:Genotype].=="Rad21KO", :]

CSV.write("Il12b_sumary_Rad21KO.csv", summary1)

summary1 = mean_burst_size_and_burst_fraction(REPS...)
summary1 = summary1[summary1[!,:Genotype].=="WT", :]

CSV.write("Il12b_sumary_WT.csv", summary1)

