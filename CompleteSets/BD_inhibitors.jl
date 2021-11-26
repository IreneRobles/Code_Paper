AnalysisName= "Il12bHSS1"

folder = "../TSS_quantification/"*AnalysisName*"/"
gene = "Il12b_intron_BD"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
e3 = get_data(folder, gene, 3)
e4 = get_data(folder, gene, 4)
e5 = get_data(folder, gene, "5_37")
e6 = get_data(folder, gene, "6_60")
e7 = get_data(folder, gene, "7_60")
exps = join_in_all_common_columns(e1,e2, e3, e4, e5, e6, e7)
CSV.write("GeneData/"*gene*".csv", exps)


folder = "../TSS_quantification/"*AnalysisName*"/"
gene = "HSS1_BD"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
e3 = get_data(folder, gene, 3)
e4 = get_data(folder, gene, 4)
e5 = get_data(folder, gene, "5_30")
e6 = get_data(folder, gene, "6_70")
e7 = get_data(folder, gene, "7_47")
exps = join_in_all_common_columns(e1,e2, e3, e4, e5, e6, e7)
CSV.write("GeneData/"*gene*".csv", exps)

folder = "../TSS_quantification/"*AnalysisName*"/"
gene = "Egr2_enh_BD"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
e3 = get_data(folder, gene, 3)
e4 = get_data(folder, gene, 4)
e5 = get_data(folder, gene, "5_37")
e6 = get_data(folder, gene, "6_60")
e7 = get_data(folder, gene, "7_60")
exps = join_in_all_common_columns(e1,e2, e3, e4, e5, e6, e7)
CSV.write("GeneData/"*gene*".csv", exps)


folder = "../TSS_quantification/"*AnalysisName*"/"
gene = "Egr2_intron_BD"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
e3 = get_data(folder, gene, 3)
e4 = get_data(folder, gene, 4)
e5 = get_data(folder, gene, "5_30")
e6 = get_data(folder, gene, "6_70")
e7 = get_data(folder, gene, "7_47")
exps = join_in_all_common_columns(e1,e2, e3, e4, e5, e6, e7)
CSV.write("GeneData/"*gene*".csv", exps)

