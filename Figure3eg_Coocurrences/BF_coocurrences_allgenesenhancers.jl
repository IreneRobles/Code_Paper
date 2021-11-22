suff = "1um_noz"
genefolder = "../CompleteSets/linkedlocus"
nascent = "Egr2_intron"
ehn = "Egr2_enh"
egr2 = DataFrame(CSV.File(normpath(genefolder, nascent*"__"*ehn*"__linked"*suff*".csv")))

nascent = "Prdm1_intron"
ehn = "Prdm1_enh"
prdm1 = DataFrame(CSV.File(normpath(genefolder, nascent*"__"*ehn*"__linked"*suff*".csv")))
prdm1 = prdm1[[in(ii, [0, 60]) for ii in prdm1[!,:Timepoint]], :]

nascent = "Il12b_intron"
ehn = "HSS1"
il12b = DataFrame(CSV.File(normpath(genefolder, nascent*"__"*ehn*"__linked"*suff*".csv")))

nascent = "Ifnb1forL2"
ehn = "L2"
ifnb1 = DataFrame(CSV.File(normpath(genefolder, nascent*"__"*ehn*"__linked"*suff*".csv")))

nascent = "Peli1_intron"
ehn = "Enh"
peli1 = DataFrame(CSV.File(normpath(genefolder, nascent*"__"*ehn*"__linked"*suff*".csv")))

nascent = "Sertad2_intron"
ehn = "Enh"
sertad2 = DataFrame(CSV.File(normpath(genefolder, nascent*"__"*ehn*"__linked"*suff*".csv")))

probe1= "Il12b_intron"
probe2 = "HSS1"
tb = il12b
cooc_il12b = DataFrame(CSV.File(normpath(genefolder, probe1*"__"*probe2*"__coocurrences"*suff*".csv")))

probe1= "Prdm1_intron"
probe2 = "Prdm1_enh"
tb = prdm1
cooc_prdm1 = DataFrame(CSV.File(normpath(genefolder, probe1*"__"*probe2*"__coocurrences"*suff*".csv")))
cooc_prdm1[!,:Timepoint] = [split(ii, "_")[2] for ii in cooc_prdm1[!,:Sample_Rep]]
cooc_prdm1 = cooc_prdm1[[in(ii, ["0", "60"]) for ii in cooc_prdm1[!,:Timepoint]], :]

probe1= "Egr2_intron"
probe2 = "Egr2_enh"
tb = egr2
cooc_egr2 = DataFrame(CSV.File(normpath(genefolder, probe1*"__"*probe2*"__coocurrences"*suff*".csv")))

probe1 = "Ifnb1forL2"
probe2 = "L2"
tb = ifnb1
cooc_ifnb1 = DataFrame(CSV.File(normpath(genefolder, probe1*"__"*probe2*"__coocurrences"*suff*".csv")))

probe1 = "Peli1_intron"
probe2 = "Enh"
tb = peli1
cooc_peli = DataFrame(CSV.File(normpath(genefolder, probe1*"__"*probe2*"__coocurrences"*suff*".csv")))

probe1 = "Sertad2_intron"
probe2 = "Enh"
tb = sertad2
cooc_sertad = DataFrame(CSV.File(normpath(genefolder, probe1*"__"*probe2*"__coocurrences"*suff*".csv")))

ts = []

probe1= "Egr2_intron"
probe2 = "Egr2_enh"
t = cooc_egr2

t2 = t
t2[!,"BF_Enhancer"] = t2[!,"BF_"*probe2]
t2[!,"BF_Gene"] = t2[!,"BF_"*probe1]
t2[!,"Gene"] = ["Egr2" for ii in 1:nrow(t2)]

push!(ts, t2)

probe1= "Il12b_intron"
probe2 = "HSS1"
t = cooc_il12b

t2 = t
t2[!,"BF_Enhancer"] = t2[!,"BF_"*probe2]
t2[!,"BF_Gene"] = t2[!,"BF_"*probe1]
t2[!,"Gene"] = ["Il12b" for ii in 1:nrow(t2)]

push!(ts, t2)

probe1= "Prdm1_intron"
probe2 = "Prdm1_enh"
t = cooc_prdm1

t2 = t
t2[!,"BF_Enhancer"] = t2[!,"BF_"*probe2]
t2[!,"BF_Gene"] = t2[!,"BF_"*probe1]
t2[!,"Gene"] = ["Prdm1" for ii in 1:nrow(t2)]

push!(ts, t2)

probe1= "Ifnb1forL2"
probe2 = "L2"
t = cooc_ifnb1

t2 = t
t2[!,"BF_Enhancer"] = t2[!,"BF_"*probe2]
t2[!,"BF_Gene"] = t2[!,"BF_"*probe1]
t2[!,"Gene"] = ["Ifnb1" for ii in 1:nrow(t2)]

push!(ts, t2)

probe1= "Peli1_intron"
probe2 = "Enh"
t = cooc_peli

t2 = t
t2[!,"BF_Enhancer"] = t2[!,"BF_"*probe2]
t2[!,"BF_Gene"] = t2[!,"BF_"*probe1]
t2[!,"Gene"] = ["Peli1" for ii in 1:nrow(t2)]

push!(ts, t2)


tab = join_in_all_common_columns(ts...)
CSV.write("BF_coocurrences_allgenesenhancers.csv", tab)