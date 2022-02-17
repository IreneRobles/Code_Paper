
nopoints = ["PAN-LPS", "BD2-5u-LPS"]

AnalysisName= "Il12bHSS1"

folder = "../TSS_quantification/"*AnalysisName*"/"
gene = "Il12b_intron_BD"
e1 = get_data(folder, gene, 1)
#very few cells
e1 = e1[e1[!,:Sample].!="Rad21KO_DMSO-LPS", :]
e2 = get_data(folder, gene, 2)
e5 = get_data2(folder, gene, "5_37")
e6 = get_data(folder, gene, "6_60")
e7 = get_data(folder, gene, "7_60")
exps = join_in_all_common_columns(e1,e2, e5, e6, e7)
exps = exps[[!in(ii, nopoints) for ii in exps[!, :Timepoint]], :]
#high background
bool1 = exps[!,:Sample].=="Rad21KO_DMSO-LPS"
bool2 = exps[!,:Rep].==8
exps = exps[.!(bool1.*bool2), :]
CSV.write("GeneData/"*gene*".csv", exps)


folder = "../TSS_quantification/"*AnalysisName*"/"
gene = "HSS1_BD"
e1 = get_data(folder, gene, 1)
#very few cells
e1 = e1[e1[!,:Sample].!="Rad21KO_DMSO-LPS", :]
e2 = get_data(folder, gene, 2)
e5 = get_data2(folder, gene, "5_30")
e6 = get_data(folder, gene, "6_70")
e7 = get_data(folder, gene, "7_47")
exps = join_in_all_common_columns(e1,e2, e5, e6, e7)
exps = exps[[!in(ii, nopoints) for ii in exps[!, :Timepoint]], :]
#high background
bool1 = exps[!,:Sample].=="Rad21KO_DMSO-LPS"
bool2 = exps[!,:Rep].==8
exps = exps[.!(bool1.*bool2), :]
exps = exps[[!in(ii, nopoints) for ii in exps[!, :Timepoint]], :]

CSV.write("GeneData/"*gene*".csv", exps)

AnalysisName= "Egr2"

folder = "../TSS_quantification/"*AnalysisName*"/"
gene = "Egr2_enh_BD"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
e3 = get_data(folder, gene, 3)
e4 = get_data(folder, gene, 4)

exps = join_in_all_common_columns(e1,e2, e4)
exps = exps[[!in(ii, nopoints) for ii in exps[!, :Timepoint]], :]

CSV.write("GeneData/"*gene*".csv", exps)

folder = "../TSS_quantification/"*AnalysisName*"/"
gene = "Egr2_intron_BD"
e1 = get_data(folder, gene, 1)
e2 = get_data(folder, gene, 2)
e3 = get_data(folder, gene, 3)
e4 = get_data(folder, gene, 4)

exps = join_in_all_common_columns(e1,e2, e4)
exps = exps[[!in(ii, nopoints) for ii in exps[!, :Timepoint]], :]

CSV.write("GeneData/"*gene*".csv", exps)

