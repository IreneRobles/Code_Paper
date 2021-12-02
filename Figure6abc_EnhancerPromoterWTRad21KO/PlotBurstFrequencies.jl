r= 3
c = 5
f = plot_bf
ii = 0
enhlabel = "Enh. burst freq."
figure(figsize = (12, 6))

subplot(r, c, ii+=1)
gene = "HSS1"
hs = [0.12, 0.20]
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
tests = [do_mantelhaen(genefreq, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq, "WT_90", "Rad21KO_90")]
add_tests3(tests, hs, u = maximum(hs)/20)
ylabel(enhlabel)


subplot(r, c, ii+=1)
gene = "Prdm1_enh"
hs = [0.12, 0.20]
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
tests = [do_mantelhaen(genefreq, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq, "WT_60", "Rad21KO_60")]
add_tests3(tests, hs, u = maximum(hs)/20)
ylabel(enhlabel)

subplot(r, c, ii+=1)
gene = "Enh"
hs = [0.12, 0.40]
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
tests = [do_mantelhaen(genefreq, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq, "WT_120", "Rad21KO_120")]
add_tests3(tests, hs, u = maximum(hs)/20)
ylabel(enhlabel)
title("Peli1_enh")



subplot(r, c, ii+=1)
gene = "Egr2_enh"
hs = [0.25, 0.35]
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
tests = [do_mantelhaen(genefreq, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq, "WT_60", "Rad21KO_60")]
add_tests3(tests, hs, u = maximum(hs)/20)
ylabel(enhlabel)


subplot(r, c, ii+=1)
gene = "L2"
hs = [0.03, 0.04]
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
tests = [do_mantelhaen(genefreq, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq, "WT_90", "Rad21KO_90")]
add_tests3(tests, hs, u = maximum(hs)/20)
ylabel(enhlabel)


subplot(r, c, ii+=1)
gene = "Il12b"
hs = [0.05, 0.13]
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
tests = [do_mantelhaen(genefreq, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq, "WT_90", "Rad21KO_90")]
add_tests3(tests, hs, u = maximum(hs)/20)
ylabel(enhlabel)


subplot(r, c, ii+=1)
gene = "Prdm1"
hs = [0.05, 0.06]
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
tests = [do_mantelhaen(genefreq, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq, "WT_60", "Rad21KO_60")]
add_tests3(tests, hs, u = maximum(hs)/20)
ylabel(enhlabel)

subplot(r, c, ii+=1)
gene = "Peli1"
hs = [0.45, 0.95]
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
tests = [do_mantelhaen(genefreq, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq, "WT_120", "Rad21KO_120")]
add_tests3(tests, hs, u = maximum(hs)/20)
ylabel(enhlabel)
title("Peli1_enh")



subplot(r, c, ii+=1)
gene = "Egr2"
hs = [0.15, 0.45]
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
tests = [do_mantelhaen(genefreq, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq, "WT_60", "Rad21KO_60")]
add_tests3(tests, hs, u = maximum(hs)/20)
ylabel(enhlabel)


subplot(r, c, ii+=1)
gene = "Ifnb1"
hs = [0.03, 0.156]
f(gene, scale = "min")
genefreq = BFs[BFs[!,:Gene].==gene, :]
tests = [do_mantelhaen(genefreq, "WT_0", "Rad21KO_0"),
    do_mantelhaen(genefreq, "WT_90", "Rad21KO_90")]
add_tests3(tests, hs, u = maximum(hs)/20)
ylabel(enhlabel)



