ehnname = "HSS1_timecourse"
intronname = "Il12b_intron_timecourse"
AnalysisName="Il12bHSS1/"

folder = "../TSS_quantification/"*AnalysisName*"/"

exp1_enh = get_data(folder, ehnname, 1;)
exp2_enh  = get_data(folder, ehnname, 2;)
exp1_gene = get_data(folder, intronname, 1;)
exp2_gene  = get_data(folder, intronname, 2;);

all_ehn = join_in_all_common_columns(exp1_enh, exp2_enh)
CSV.write("GeneData/"*"$ehnname"*".csv", all_ehn)
all_gene = join_in_all_common_columns(exp1_gene, exp2_gene)
CSV.write("GeneData/"*"$intronname"*".csv", all_gene)