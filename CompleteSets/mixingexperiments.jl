folder = "Ifit1/"
folder = "../TSS_quantification/"*folder*"/"
gene = "Ifit1"
e1 = get_data(folder, gene, 1, addrep = 1)
e2 = get_data(folder, gene, 2, addrep = 2)
e3 = get_data(folder, gene, 3, addrep = 3)
exps = vcat(e1,e2, e3)
wt = exps[exps[!,:Genotype].=="WT", :]
rad = exps[exps[!,:Genotype].=="Rad21KO", :]
exps = vcat(wt, rad)
CSV.write("GeneData/"*gene*".csv", exps)

folder = "Cxcl10/"
folder = "../TSS_quantification/"*folder*"/"

gene = "Cxcl10"
e1 = get_data(folder, gene, 1, addrep = 1)
e2 = get_data(folder, gene, 2, addrep = 2)
e3 = get_data(folder, gene, 3, addrep = 3)
e4 = get_data(folder, gene, 4, addrep = 4)
exps = vcat(e1,e2, e3, e4)
wt = exps[exps[!,:Genotype].=="WT", :]
rad = exps[exps[!,:Genotype].=="Rad21KO", :]
exps = vcat(wt, rad)
CSV.write("GeneData/"*gene*".csv", exps)