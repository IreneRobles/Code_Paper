
AnalysisName= "Il12bHSS1/"

folder = "../TSS_quantification/"*AnalysisName*"/"

gene = "Il12b_intron_IFNG"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)

CSV.write("GeneData/"*gene*".csv", exps)


AnalysisName= "Il12bHSS1/"
folder = "../TSS_quantification/"*AnalysisName*"/"

gene = "HSS1_IFNG"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)

CSV.write("GeneData/"*gene*".csv", exps)



AnalysisName= "Egr2/"

folder = "../TSS_quantification/"*AnalysisName*"/"

gene = "Egr2_intron_IFNG"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)

CSV.write("GeneData/"*gene*".csv", exps)


AnalysisName= "Egr2/"
folder = "../TSS_quantification/"*AnalysisName*"/"

gene = "Egr2_enh_IFNG"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
exps = vcat(e1,e2)

CSV.write("GeneData/"*gene*".csv", exps)