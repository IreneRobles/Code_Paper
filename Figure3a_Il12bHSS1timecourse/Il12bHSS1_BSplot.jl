
limity = 28
u = limity/20
subplot(3,2,3)
genesize_gene = plot_burstsizes(all_gene, GENE ="$intronname", limit = 1, palette = ["darkgray"], order = nothing)
filename = "$AnalysisName/"*"$intronname"*"_BS_Test.csv"
R"""
tb <- $genesize_gene
tb$Timepoint <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Timepoint +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Timepoint
write.csv(a, $filename)
"""
t = CSV.read(filename, DataFrames.DataFrame)


s1 = "WT_0"
s2 = "WT_60"
h = 17.5
cords = [0,1]
add_test(t, s1, s2, h, u, cords)


s1 = "WT_60"
s2 = "WT_90"
p = 1
h = 20.0
cords = [1,2]
add_test(t, s1, s2, h, u, cords)

s1 = "WT_90"
s2 = "WT_120"
h = 17.0
cords = [2,3]
add_test(t, s1, s2, h, u, cords)

s1 = "WT_60"
s2 = "WT_120"
h = 23.0
cords = [1,3]
add_test(t, s1, s2, h, u, cords)


ylim(0,limity);legend_removal();pretty_axes2();line075black();squareplot();xticks(rotation=0), title("Il12b")
ylabel("Gene burst size")



subplot(3,2,4)

limity = 8.0
u = limity/20

genesize_gene = plot_burstsizes(all_ehn, GENE ="$ehnname", limit = 1, palette = ["darkgray"], order = nothing)
filename = "$AnalysisName/"*"$intronname"*"_BS_Test.csv"
R"""
tb <- $genesize_gene
tb$Timepoint <- as.factor(tb$Sample)
tb$Rep <- as.factor(tb$Rep); 
tb$TSS1_r2 <- as.numeric(tb$TSS1_r2); 
aov.s = aov(TSS1_r2 ~ Timepoint +Rep ,data=tb)  #do the analysis of variance
test <- TukeyHSD(aov.s)
a = test$Timepoint
write.csv(a, $filename)
"""
t = CSV.read(filename, DataFrames.DataFrame)

s1 = "WT_0"
s2 = "WT_60"
h = 5.0
cords = [0,1]
add_test(t, s1, s2, h, u, cords)


s1 = "WT_60"
s2 = "WT_120"
p = 1
h = 7.0
cords = [1,3]
add_test(t, s1, s2, h, u, cords)

s1 = "WT_60"
s2 = "WT_90"
p = 1
h = 6.0
cords = [1,2]
add_test(t, s1, s2, h, u, cords)

s1 = "WT_90"
s2 = "WT_120"
h = 5.0
cords = [2,3]
add_test(t, s1, s2, h, u, cords)


ylim(0,limity);legend_removal();pretty_axes2();line075black();squareplot();xticks(rotation=0); title("HSS1")
ylabel("Enh burst size")
plt.tight_layout()