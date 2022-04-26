


function plotgeneenhancerlog2bf(gene; Genotype = "All")
a = tb
    print("gene")

if Genotype != "All"
a = a[a[!,:Genotype].==Genotype, :]
end

figure(figsize = [3,3])
a[!,"log2_BF_Enhancer"] = log2.(a[!,"BF_Enhancer"])
a[!,"log2_BF_Gene"] = log2.(a[!,"BF_Gene"]) 
    

a = a[.! isinf.(a[!,"log2_BF_Gene"]), :]
a = a[.! isinf.(a[!,"log2_BF_Enhancer"]), :]    

x = a[!,"log2_BF_Enhancer"]
y = a[!,"log2_BF_Gene"]


pd = Pandas.DataFrame(a)
pd["LPS"] = pd["Timepoint"] .> 0

    
    
if Genotype != "All"
py"""
import seaborn as Seaborn
Seaborn.scatterplot(data = $pd, x = "log2_BF_Enhancer",  y = "log2_BF_Gene", hue = "Gene", s = 100, style = "LPS")
Seaborn.regplot(data = $pd, x = "log2_BF_Enhancer",  y = "log2_BF_Gene", scatter = 0, color = "grey", ci =  95 )
"""
else
        py"""
import seaborn as Seaborn
Seaborn.scatterplot(data = $pd, x = "log2_BF_Enhancer",  y = "log2_BF_Gene", hue = "Genotype",palette = ["black", "red"], s = 100, style = "Gene")
b = $pd
wt = b[b["Genotype"]=="WT"]
    
Seaborn.regplot(data = wt, x = "log2_BF_Enhancer",  y = "log2_BF_Gene", scatter = 0, color = "black", ci =  95 )

rad = b[b["Genotype"]!="WT"]
    
Seaborn.regplot(data = rad, x = "log2_BF_Enhancer",  y = "log2_BF_Gene", scatter = 0, color = "red", ci =  95 )

"""
end


TEST = R"""cor.test($x,$y)"""
corr= round(TEST[4][1], digits = 3)
pval= round(TEST[3][1], sigdigits = 3)

    annotate("r = $corr \np = $pval", xy = (-10, -2))
title("$gene")
ylim(-12 , 0.5)
    
    xlim(-10, 0.5)

    ax = gca()
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable="box")



ylabel("log2(Gene \n Burst Frequency)")

xlabel("log2(Enhancer \n Burst Frequency)")
    
    legend_out_of_plot()
    pretty_axes2()

    
    
    py"""
import matplotlib as mpl

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
mpl.rcParams.update(new_rc_params)
import matplotlib.pyplot as plt
plt.savefig($gene+"_WT_BFenhancerBFgene.svg")
"""
    TEST
end

function plotgeneenhancerbf(gene, tb; xy = (0.01, -0.01), onlyLPS = false)
a = tb
println(gene)
a = a[a[!,"Gene"].== gene, :]

if onlyLPS
        a = a[a[!,:Timepoint].!= 0, :]
    end
    
if gene == "Prdm1_intron"
        a = a[a[!,:Timepoint].!= "30", :]
        a = a[a[!,:Timepoint].!= "90", :]
        
    end
    

pd = Pandas.DataFrame(a)
pd["LPS"] = pd["Timepoint"] .> 0
    


    

py"""
import seaborn as Seaborn
Seaborn.scatterplot(data = $pd, x = "BF_Enhancer",  y = "BF_Gene", hue = "Genotype", style = "LPS",palette = ["black", "red"], s = 100)
#Seaborn.scatterplot(data = $pd, x = "BF_Enhancer",  y = "BF_Gene", hue = "Rep", s = 100)

b = $pd
wt = b[b["Genotype"]=="WT"]
    
Seaborn.regplot(data = wt, x = "BF_Enhancer",  y = "BF_Gene", scatter = 0, color = "black", ci =  95 )

rad = b[b["Genotype"]!="WT"]
    
Seaborn.regplot(data = rad, x = "BF_Enhancer",  y = "BF_Gene", scatter = 0, color = "red", ci =  95 )

"""

b = a[a[!,:Genotype].== "WT", :]
x = b[!,"BF_Enhancer"]
y = b[!,"BF_Gene"]

TEST = R"""cor.test($x,$y)"""
corr= round(TEST[4][1], digits = 3)
pval= round(TEST[3][1], sigdigits = 3)

    print("WT\nr = $corr\np = $pval")
    
b = a[a[!,:Genotype].== "Rad21KO", :]
x = b[!,"BF_Enhancer"]
y = b[!,"BF_Gene"]

TEST = R"""cor.test($x,$y)"""
corr= round(TEST[4][1], digits = 3)
pval= round(TEST[3][1], sigdigits = 3)

    print("\nRad21KO\nr = $corr\np = $pval")
    
b = a
x = b[!,"BF_Enhancer"]
y = b[!,"BF_Gene"]

TEST = R"""cor.test($x,$y)"""
corr= round(TEST[4][1], digits = 3)
pval= round(TEST[3][1], sigdigits = 3)

    print("\nAll\nr = $corr\np = $pval")

title("$gene")

    ax = gca()
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable="box")



ylabel("Gene Burst \n Frequency")

xlabel("Enhancer Burst \n Frequency")
    
    legend_out_of_plot()
    pretty_axes2()


    R"""
tb <- $a
testlm <- lm(BF_Gene ~ BF_Enhancer*Genotype, data = tb)
print(summary(aov(testlm)))
"""
    
end






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
figure(figsize = (15, 8))

subplot(r, c, 1)
gene = "HSS1"
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
sort!(genefreq,[:Rep, :Sample] )
genefreq_enh = genefreq
tests = [do_mantelhaen(genefreq_enh, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_enh, "WT_90", "Rad21KO_90")]
ylim(0, 0.4)
add_tests3(tests, [0.1, 0.2], u = 0.02)

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
add_tests3(tests, [0.2, 0.3].*n, u = 0.02*n)

plt.tight_layout()

subplot(r, c, 12)

plot_BFGene_Enhancer(genefreq_gene, genefreq_enh)

subplot(r, c, 3)
n = 1
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

n = 2
gene = "Peli1"
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
sort!(genefreq,[:Rep, :Sample] )
genefreq_gene = genefreq[genefreq[!,:Rep].>3, :]


tests = [do_mantelhaen(genefreq_gene, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_gene, "WT_120", "Rad21KO_120")]
ylim(0, 0.4*n)
add_tests3(tests, [0.22, 0.5].*n, u = 0.02*n)


subplot(r, c, 13)
plot_BFGene_Enhancer(genefreq_gene, genefreq_enh)


subplot(r, c, 4)
n = 2
gene = "Egr2_enh"
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]

sort!(genefreq,[:Rep, :Sample] )
genefreq_enh = genefreq
tests = [do_mantelhaen(genefreq_enh, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq_enh, "WT_60", "Rad21KO_60")]


ylim(0, 0.4*n)
add_tests3(tests, [0.15, 0.2].*n, u = 0.02*n)

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
n = 0.2
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




#subplot(r, c, 15)

#plot_BFGene_Enhancer(genefreq_gene, genefreq_enh)


plt.tight_layout()

