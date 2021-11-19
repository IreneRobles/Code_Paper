limitbs = 1
BF = Dict()
bf = calculate_bf_by_rep(get_genedata("Enh"), "TSS", "Enh", limit = limitbs)
BF["Enh"] = bf
bf = calculate_bf_by_rep(get_genedata("L2"), "TSS", "L2", limit = limitbs)
BF["L2"] = bf
bf = calculate_bf_by_rep(get_genedata("Egr2_enh"), "TSS", "Egr2_enh", limit = limitbs)
BF["Egr2_enh"] = bf
bf = calculate_bf_by_rep(get_genedata("Prdm1_enh"), "TSS", "Prdm1_enh", limit = limitbs)
bf = bf[bf[!,:Timepoint].!= "30", :]
bf = bf[bf[!,:Timepoint].!= "90", :]
BF["Prdm1_enh"] = bf
bf = calculate_bf_by_rep(get_genedata("HSS1"), "TSS", "HSS1", limit = 0)
BF["HSS1"] = bf
BFs = join_in_all_common_columns([BF[key] for key in keys(BF)]...)
CSV.write("BurstFrequencies_Enhancers.csv", BFs)