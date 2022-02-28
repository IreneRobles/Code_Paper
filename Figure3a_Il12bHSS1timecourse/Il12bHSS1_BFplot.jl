exps = all_gene
exp1 = exps[exps[!,:Rep].==1, :]
exp2 = exps[exps[!,:Rep].==2, :];
exp3 = exps[exps[!,:Rep].==3, :]
exp4 = exps[exps[!,:Rep].==4, :];

ylimit = 0.2
u = ylimit/20

BFs_gene = calculatesumary(exp1, exp2, exp3, exp4, limit = 1)
BFs_g = BFs_gene
TEST = BF_test(BFs_gene, "$intronname", :N_TSS, [
        ["WT_0", "WT_60"],
        ["WT_60", "WT_90"],
        ["WT_90", "WT_120"]
        ])

CSV.write("$AnalysisName/"*"$intronname"*"_BFs.csv", BFs_gene)
#CSV.write("$AnalysisName/"*"$intronname"*"_BFs_test.csv", TEST)
figure(figsize = (6, 8))
subplot(3,2,1);title("Il12b")
    Seaborn.stripplot(data = Pandas.DataFrame(BFs_gene), x = "Timepoint", y = "BF", hue = "Rep", palette = "Greys")
Seaborn.boxplot(data = Pandas.DataFrame(BFs_gene), y = "BF", x = "Timepoint", hue = "Genotype", dodge = true, palette = ["darkgray"], showfliers = false)

h1,h2,h3, = 0.05, 0.14, 0.17

ii = 0
h = h1
xy = [ii, ii+1]; xy2 = [Statistics.mean(xy) , h+u]
plt.plot(xy, [h, h], c = "black")
pval = TEST[ii+=1, :padj]
plt.annotate(NoLongerProblems.transform_pvalue_in_stars(pval), xy = xy2, ha= "center")

h = h2
xy = [ii, ii+1]; xy2 = [Statistics.mean(xy) , h+u]
plt.plot(xy, [h, h], c = "black")
pval = TEST[ii+=1, :padj]
plt.annotate(NoLongerProblems.transform_pvalue_in_stars(pval), xy = xy2, ha= "center")

h = h3
xy = [ii, ii+1]; xy2 = [Statistics.mean(xy) , h+u]
plt.plot(xy, [h, h], c = "black")
pval = TEST[ii+=1, :padj]
plt.annotate(NoLongerProblems.transform_pvalue_in_stars(pval), xy = xy2, ha= "center")
ylim(0,ylimit)
xlabel("Time after LPS (min)")
ylabel("Gene burst freq.");legend_removal();pretty_axes2();line075black();squareplot();xticks(rotation=0)

exps = all_ehn
exp1 = exps[exps[!,:Rep].==1, :]
exp2 = exps[exps[!,:Rep].==2, :];
exp3 = exps[exps[!,:Rep].==3, :]
exp4 = exps[exps[!,:Rep].==4, :];

ylimit = 0.5
u = ylimit/20

BFs_gene = calculatesumary(exp1, exp2, exp3, exp4, limit = 0)
BFs_e = BFs_gene


TEST = BF_test(BFs_gene, "$ehnname", :N_TSS, [
        ["WT_0", "WT_60"],
        ["WT_60", "WT_90"],
        ["WT_90", "WT_120"]
        ])

CSV.write("$AnalysisName/"*"$ehnname"*"_BFs.csv", BFs_gene)
#CSV.write("$AnalysisName/"*"$ehnname"*"_BFs_test.csv", TEST)

subplot(3,2,2);title("HSS1")
    Seaborn.stripplot(data = Pandas.DataFrame(BFs_gene), x = "Timepoint", y = "BF", hue = "Rep", palette = "Greys")
Seaborn.boxplot(data = Pandas.DataFrame(BFs_gene), y = "BF", x = "Timepoint", hue = "Genotype", dodge = true, palette = ["darkgray"])

h1,h2,h3, = 0.45, 0.4, 0.3

ii = 0
h = h1
xy = [ii, ii+1]; xy2 = [Statistics.mean(xy) , h+u]
plt.plot(xy, [h, h], c = "black")
pval = TEST[ii+=1, :padj]
plt.annotate(NoLongerProblems.transform_pvalue_in_stars(pval), xy = xy2, ha= "center")

h = h2
xy = [ii, ii+1]; xy2 = [Statistics.mean(xy) , h+u]
plt.plot(xy, [h, h], c = "black")
pval = TEST[ii+=1, :padj]
plt.annotate(NoLongerProblems.transform_pvalue_in_stars(pval), xy = xy2, ha= "center")

h = h3
xy = [ii, ii+1]; xy2 = [Statistics.mean(xy) , h+u]
plt.plot(xy, [h, h], c = "black")
pval = TEST[ii+=1, :padj]
plt.annotate(NoLongerProblems.transform_pvalue_in_stars(pval), xy = xy2, ha= "center")
ylim(0,ylimit)
ylabel("Enh burst freq.");legend_removal();pretty_axes2();line075black();squareplot();xticks(rotation=0)
xlabel("Time after LPS (min)")
plt.tight_layout()

