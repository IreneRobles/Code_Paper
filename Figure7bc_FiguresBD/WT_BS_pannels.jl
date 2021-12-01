subplot(2,4,1)  
bs = bs_egr2enh
gen = "WT"
titlegene = "Egr2 Enh"
bs_test = CSV.read("Egr2_enh_BD_BS_Test.csv", DataFrames.DataFrame)
n = 1
ylimit = 7*n
u = ylimit/30
plot_burstsizes_genotype(bs, gen, ylimit = ylimit, limit = hss1lim, GENE = titlegene)
findpadj(bs_test, gen*"_DMSO",gen*"_DMSO-LPS", xy = (0,1), h = 5.5*n, u = u)
findpadj(bs_test, gen*"_DMSO-LPS",gen*"_BD1-LPS", xy = (1,2), h = 5*n, u = u)
findpadj(bs_test, gen*"_DMSO-LPS",gen*"_BD2-LPS", xy = (1,3), h = 6*n, u = u)
findpadj(bs_test, gen*"_BD1-LPS",gen*"_BD2-LPS", xy = (2,3), h = 5.25*n, u = u)
ylabel("Enh. burst size")
    
subplot(2,4,2)  
bs = bs_egr2
gen = "WT"
titlegene = "Egr2"
bs_test = CSV.read("Egr2_intron_BD_BS_Test.csv", DataFrames.DataFrame)
n = 1.5
ylimit = 7*n
u = ylimit/30
plot_burstsizes_genotype(bs, gen, ylimit = ylimit, limit = 1, GENE = titlegene)
findpadj(bs_test, gen*"_DMSO",gen*"_DMSO-LPS", xy = (0,1), h = 5.5*n, u = u)
findpadj(bs_test, gen*"_DMSO-LPS",gen*"_BD1-LPS", xy = (1,2), h = 5*n, u = u)
findpadj(bs_test, gen*"_DMSO-LPS",gen*"_BD2-LPS", xy = (1,3), h = 6*n, u = u)
findpadj(bs_test, gen*"_BD1-LPS",gen*"_BD2-LPS", xy = (2,3), h = 5.25*n, u = u)
ylabel("Gene burst size")


subplot(2,4,3)  
bs = bs_hss1
gen = "WT"
titlegene = "Il12b HSS1 Enh"
bs_test = CSV.read("HSS1_BD_BS_Test.csv", DataFrames.DataFrame)
n = 0.5
ylimit = 7*n
u = ylimit/30
plot_burstsizes_genotype(bs, gen, ylimit = ylimit, limit = hss1lim, GENE = titlegene)
findpadj(bs_test, gen*"_DMSO",gen*"_DMSO-LPS", xy = (0,1), h = 5.5*n, u = u)
findpadj(bs_test, gen*"_DMSO-LPS",gen*"_BD1-LPS", xy = (1,2), h = 5*n, u = u)
findpadj(bs_test, gen*"_DMSO-LPS",gen*"_BD2-LPS", xy = (1,3), h = 6*n, u = u)
findpadj(bs_test, gen*"_BD1-LPS",gen*"_BD2-LPS", xy = (2,3), h = 5.25*n, u = u)
ylabel("Enh. burst size")

    
subplot(2,4,4)  
bs = bs_il12b
gen = "WT"
titlegene = "Il12b"
bs_test = CSV.read("Il12b_intron_BD_BS_Test.csv", DataFrames.DataFrame)
n = 3.3
ylimit = 7*n
u = ylimit/30
plot_burstsizes_genotype(bs, gen, ylimit = ylimit, limit = 1, GENE = titlegene)
findpadj(bs_test, gen*"_DMSO",gen*"_DMSO-LPS", xy = (0,1), h = 5.5*n, u = u)
findpadj(bs_test, gen*"_DMSO-LPS",gen*"_BD1-LPS", xy = (1,2), h = 5*n, u = u)
findpadj(bs_test, gen*"_DMSO-LPS",gen*"_BD2-LPS", xy = (1,3), h = 6*n, u = u)
findpadj(bs_test, gen*"_BD1-LPS",gen*"_BD2-LPS", xy = (2,3), h = 5.25*n, u = u)
ylabel("Gene burst size")