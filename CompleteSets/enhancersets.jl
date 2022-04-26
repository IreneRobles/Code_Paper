folder = "EhnPeliSertad/"
gene = "Enh"
folder = "../TSS_quantification/"*folder*"/"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)
CSV.write("GeneData/"*gene*".csv", exps)
gene = "Peli1_intron"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)
CSV.write("GeneData/"*gene*".csv", exps)
gene = "Sertad2_intron"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)
CSV.write("GeneData/"*gene*".csv", exps)

ehnname = "HSS1"
intronname = "Il12b_intron"
AnalysisName="Il12bHSS1/"

folder = "../TSS_quantification/"*AnalysisName*"/"

exp1_enh = get_data(folder, ehnname, "1_extra";)
exp2_enh  = get_data(folder, ehnname, 2;)
exp1_gene = get_data(folder, intronname, "1_extra";)
exp2_gene  = get_data(folder, intronname, 2;);

all_ehn = join_in_all_common_columns(exp1_enh, exp2_enh)
CSV.write("GeneData/"*"$ehnname"*".csv", all_ehn)
all_gene = join_in_all_common_columns(exp1_gene, exp2_gene)
CSV.write("GeneData/"*"$intronname"*".csv", all_gene)

folder = "Egr2"
gene = "Egr2_intron"
folder = "../TSS_quantification/"*folder*"/"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
#e2[!,:Rep] = [2 for ii in 1:nrow(e2)]

exps = vcat(e1,e2)
CSV.write("GeneData/"*gene*".csv", exps)
gene = "Egr2_enh"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
#e2[!,:Rep] = [2 for ii in 1:nrow(e2)]

exps = vcat(e1,e2)
CSV.write("GeneData/"*gene*".csv", exps)

folder = "Ifnb1L2"
gene = "Ifnb1forL2"
folder = "../TSS_quantification/"*folder*"/"

e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)

CSV.write("GeneData/"*gene*".csv", exps)

folder = "Ifnb1L2"
gene = "L2"
folder = "../TSS_quantification/"*folder*"/"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)
CSV.write("GeneData/"*gene*".csv", exps)


folder = "Prdm1/"
folder = "../TSS_quantification/"*folder*"/"

gene = "Prdm1_intron"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)
CSV.write("GeneData/"*gene*".csv", exps)
gene = "Prdm1_enh"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)
CSV.write("GeneData/"*gene*".csv", exps)