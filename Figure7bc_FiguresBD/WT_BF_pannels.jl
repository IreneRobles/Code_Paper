subplot(2,4,5)
n = 1
plot_BF_Genotype(BFs_all_gene_egr2enh,TEST_gene_egr2enh;
    hs = [0.4, 0.55, 0.63, 0.50].*n,
    ylimit = 0.7*n, 
    Genotype = "WT",
    colorstripplot = "gray",
    colorboxplot = "darkgray",
    ylabel_ = "Enh. burst frequency")

subplot(2,4,6)
n = 1
plot_BF_Genotype(BFs_all_gene_egr2,TEST_gene_egr2;
    hs = [0.4, 0.55, 0.63, 0.50].*n,
    ylimit = 0.7*n, 
    Genotype = "WT",
    colorstripplot = "gray",
    colorboxplot = "darkgray",
    ylabel_ = "Gene burst frequency")

subplot(2,4,7)
n = 0.6
plot_BF_Genotype(BFs_all_gene_hss1,TEST_gene_hss1;
    hs = [0.4, 0.55, 0.63, 0.50].*n,
    ylimit = 0.7*n, 
    Genotype = "WT",
    colorstripplot = "gray",
    colorboxplot = "darkgray",
    ylabel_ = "Enh. burst frequency")

subplot(2,4,8)
n = 0.5
plot_BF_Genotype(BFs_all_gene_il,TEST_gene_il;
    hs = [0.4, 0.55, 0.63, 0.50].*n,
    ylimit = 0.7*n, 
    Genotype = "WT",
    colorstripplot = "gray",
    colorboxplot = "darkgray",
    ylabel_ = "Gene burst frequency")