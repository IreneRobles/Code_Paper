
function plot_BFGene_Enhancer(genefreq_gene, genefreq_enh)
    t = genefreq_gene
t[!,:BF_Gene] = genefreq_gene[!,:BF_TSS]
t[!,:BF_Enh] = genefreq_enh[!,:BF_TSS]

pd = Pandas.DataFrame(t)

py"""
import seaborn
import pandas

t = $pd

seaborn.scatterplot( data = t, x= "BF_Enh", y= "BF_Gene", hue = "Genotype", style = "Timepoint",palette = ["red", "black"])

w= t[t["Genotype"]=="WT"]

seaborn.regplot( data = w, x= "BF_Enh", y= "BF_Gene",scatter = 0, color = "black")
r= t[t["Genotype"]!="WT"]

seaborn.regplot( data = r, x= "BF_Enh", y= "BF_Gene",scatter = 0, color = "red")

"""

legend_removal()
ylabel("Prom. burst freq.")
xlabel("Enh. burst freq.") 
    squareplot()
pretty_axes2()
end

r= 3
c = 5
f = plot_bf
ii = 0
enhlabel = "Enh. burst freq."
figure(figsize = (12, 6))

subplot(r, c, 1)
gene = "HSS1"
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
sort!(genefreq,[:Rep, :Sample] )
genefreq_enh = genefreq
tests = [do_mantelhaen(genefreq_enh, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_enh, "WT_90", "Rad21KO_90")]
ylim(0, 0.4)
add_tests3(tests, [0.3, 0.3], u = 0.02)

subplot(r, c, 6)
n = 0.4
gene = "Il12b"
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
sort!(genefreq,[:Rep, :Sample] )
genefreq_gene = genefreq[genefreq[!,:Rep].>4, :]
tests = [do_mantelhaen(genefreq_gene, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_gene, "WT_90", "Rad21KO_90")]
ylim(0, 0.4*n)
add_tests3(tests, [0.1, 0.3].*n, u = 0.02*n)

subplot(r, c, 11)
plot_BFGene_Enhancer(genefreq_gene, genefreq_enh)

subplot(r, c, 2)
n = 1.2
gene = "Prdm1_enh"
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]

sort!(genefreq,[:Rep, :Sample] )
genefreq_enh = genefreq
tests = [do_mantelhaen(genefreq_enh, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_enh, "WT_60", "Rad21KO_60")]


ylim(0, 0.4*n)
add_tests3(tests, [0.1, 0.3].*n, u = 0.02*n)

subplot(r, c, 7)

n = 0.2
gene = "Prdm1"
f(gene, scale = "min")
genefreq_gene = BFs[BFs[!,:Gene].==gene, :]
sort!(genefreq_gene,[:Rep, :Sample] )

tests = [do_mantelhaen(genefreq_gene, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_gene, "WT_60", "Rad21KO_60")]
ylim(0, 0.4*n)
add_tests3(tests, [0.3, 0.3].*n, u = 0.02*n)

plt.tight_layout()

subplot(r, c, 12)

plot_BFGene_Enhancer(genefreq_gene, genefreq_enh)

subplot(r, c, 3)
n = 2
gene = "Enh"
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
sort!(genefreq,[:Rep, :Sample] )
genefreq_enh = genefreq
tests = [do_mantelhaen(genefreq_enh, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_enh, "WT_120", "Rad21KO_120")]
ylim(0, 0.4*n)
add_tests3(tests, [0.1, 0.3].*n, u = 0.02*n)


subplot(r, c, 8)

n = 3
gene = "Peli1"
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
sort!(genefreq,[:Rep, :Sample] )
genefreq_gene = genefreq[genefreq[!,:Rep].>3, :]


tests = [do_mantelhaen(genefreq_gene, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_gene, "WT_120", "Rad21KO_120")]
ylim(0, 0.4*n)
add_tests3(tests, [0.22, 0.35].*n, u = 0.02*n)


subplot(r, c, 13)
plot_BFGene_Enhancer(genefreq_gene, genefreq_enh)


subplot(r, c, 4)
n = 3
gene = "Egr2_enh"
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]

sort!(genefreq,[:Rep, :Sample] )
genefreq_enh = genefreq
tests = [do_mantelhaen(genefreq_enh, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_enh, "WT_60", "Rad21KO_60")]


ylim(0, 0.4*n)
add_tests3(tests, [0.3, 0.3].*n, u = 0.02*n)

subplot(r, c, 9)

n = 2.5
gene = "Egr2"
f(gene, scale = "min")
genefreq_gene = BFs[BFs[!,:Gene].==gene, :]
sort!(genefreq_gene,[:Rep, :Sample] )

tests = [do_mantelhaen(genefreq_gene, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_gene, "WT_60", "Rad21KO_60")]
ylim(0, 0.4*n)
add_tests3(tests, [0.2, 0.3].*n, u = 0.02*n)

plt.tight_layout()



subplot(r, c, 14)

plot_BFGene_Enhancer(genefreq_gene, genefreq_enh)

subplot(r, c, 5)
n = 0.3
gene = "L2"
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
sort!(genefreq,[:Rep, :Sample] )
genefreq_enh = genefreq
tests = [do_mantelhaen(genefreq_enh, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_enh, "WT_90", "Rad21KO_90")]
ylim(0, 0.4*n)
add_tests3(tests, [0.2, 0.3].*n, u = 0.02*n)


subplot(r, c, 10)

n = 0.36
gene = "Ifnb1"

f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
sort!(genefreq,[:Rep, :Sample] )
genefreq_gene = genefreq[genefreq[!,:Rep].>1, :]
tests = [do_mantelhaen(genefreq_gene, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_gene, "WT_90", "Rad21KO_90")]
ylim(0, 0.4*n)
add_tests3(tests, [0.22, 0.35].*n, u = 0.02*n)




subplot(r, c, 15)

plot_BFGene_Enhancer(genefreq_gene, genefreq_enh)






plt.tight_layout()

