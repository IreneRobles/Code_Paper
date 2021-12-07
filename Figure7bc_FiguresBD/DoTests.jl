

# Burst frequenciy summaries and tests

BFs_all_gene_il = calculatesumary_byrep(il12b, limit = 1)
CSV.write("Il12b_BD_summary.csv", BFs_all_gene_il)

TEST_gene_il = BF_test(BFs_all_gene_il, "Il12b_intron_BD", :N_TSS, [
        ["Rad21KO_DMSO", "WT_DMSO"],
        ["Rad21KO_DMSO-LPS", "WT_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "WT_BD1-LPS"],
        ["Rad21KO_BD2-LPS", "WT_BD2-LPS"],
        ["WT_DMSO", "WT_DMSO-LPS"],
        ["WT_BD1-LPS", "WT_DMSO-LPS"],
        ["WT_BD2-LPS", "WT_DMSO-LPS"],
        ["WT_BD2-LPS", "WT_BD1-LPS"],
        ["Rad21KO_DMSO", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD2-LPS", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "Rad21KO_BD2-LPS"],
        ])

BFs_all_gene_hss1 = calculatesumary_byrep(hss1, limit = hss1lim)
CSV.write("HSS1_BD_summary.csv", BFs_all_gene_hss1)

TEST_gene_hss1 = BF_test(BFs_all_gene_hss1, "HSS1_BD", :N_TSS, [
        ["Rad21KO_DMSO", "WT_DMSO"],
        ["Rad21KO_DMSO-LPS", "WT_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "WT_BD1-LPS"],
        ["Rad21KO_BD2-LPS", "WT_BD2-LPS"],
        ["WT_DMSO", "WT_DMSO-LPS"],
        ["WT_BD1-LPS", "WT_DMSO-LPS"],
        ["WT_BD2-LPS", "WT_DMSO-LPS"],
        ["WT_BD2-LPS", "WT_BD1-LPS"],
        ["Rad21KO_DMSO", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD2-LPS", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "Rad21KO_BD2-LPS"],
        ])

BFs_all_gene_egr2 = calculatesumary_byrep(egr2, limit = 1)
CSV.write("Egr2_BD_summary.csv", BFs_all_gene_egr2)

TEST_gene_egr2 = BF_test(BFs_all_gene_egr2, "Egr2_intron_BD", :N_TSS, [
        ["Rad21KO_DMSO", "WT_DMSO"],
        ["Rad21KO_DMSO-LPS", "WT_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "WT_BD1-LPS"],
        ["Rad21KO_BD2-LPS", "WT_BD2-LPS"],
        ["WT_DMSO", "WT_DMSO-LPS"],
        ["WT_BD1-LPS", "WT_DMSO-LPS"],
        ["WT_BD2-LPS", "WT_DMSO-LPS"],
        ["WT_BD2-LPS", "WT_BD1-LPS"],
        ["Rad21KO_DMSO", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD2-LPS", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "Rad21KO_BD2-LPS"],
        ])

BFs_all_gene_egr2enh = calculatesumary_byrep(egr2enh, limit = 1)
CSV.write("Egr2enh_BD_summary.csv", BFs_all_gene_egr2enh)

TEST_gene_egr2enh = BF_test(BFs_all_gene_egr2enh, "Egr2_enh_BD", :N_TSS, [
        ["Rad21KO_DMSO", "WT_DMSO"],
        ["Rad21KO_DMSO-LPS", "WT_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "WT_BD1-LPS"],
        ["Rad21KO_BD2-LPS", "WT_BD2-LPS"],
        ["WT_DMSO", "WT_DMSO-LPS"],
        ["WT_BD1-LPS", "WT_DMSO-LPS"],
        ["WT_BD2-LPS", "WT_DMSO-LPS"],
        ["WT_BD2-LPS", "WT_BD1-LPS"],
        ["Rad21KO_DMSO", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD2-LPS", "Rad21KO_DMSO-LPS"],
        ["Rad21KO_BD1-LPS", "Rad21KO_BD2-LPS"],
        ])

# Burst size tests

bs_il12b= tss_data(il12b; limit = 1)
filename = "Il12b_intron_BD"*"_BS_Test.csv"
R"""
tb <- $bs_il12b
tb$Timepoint <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Timepoint +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Timepoint
write.csv(a, $filename)
"""
bs_il12b_test = CSV.read(filename, DataFrames.DataFrame)


bs_hss1= tss_data(hss1; limit = hss1lim)
filename = "HSS1_BD"*"_BS_Test.csv"
R"""
tb <- $bs_hss1
tb$Timepoint <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Timepoint +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Timepoint
write.csv(a, $filename)
"""
bs_hss1_test = CSV.read(filename, DataFrames.DataFrame)



bs_egr2= tss_data(egr2; limit = 1)
filename = "Egr2_intron_BD"*"_BS_Test.csv"
R"""
tb <- $bs_egr2
tb$Timepoint <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Timepoint +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Timepoint
write.csv(a, $filename)
"""
bs_egr2_test = CSV.read(filename, DataFrames.DataFrame)

bs_egr2enh= tss_data(egr2enh; limit = 1)
filename = "Egr2_enh_BD"*"_BS_Test.csv"
R"""
tb <- $bs_egr2enh
tb$Timepoint <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Timepoint +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Timepoint
write.csv(a, $filename)
"""
bs_egr2enh_test = CSV.read(filename, DataFrames.DataFrame)