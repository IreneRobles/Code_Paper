limitbs = 1
BF = Dict()
bf = calculate_bf_by_rep(get_genedata("Peli1_intron"), "TSS", "Peli1_intron", limit = limitbs)
BF["Peli1_intron"] = bf
bf = calculate_bf_by_rep(get_genedata("Ifnb1forL2"), "TSS", "Ifnb1forL2", limit = limitbs)
BF["Ifnb1forL2"] = bf
bf = calculate_bf_by_rep(get_genedata("Egr2_intron"), "TSS", "Egr2_intron", limit = limitbs)
BF["Egr2_intron"] = bf
bf = calculate_bf_by_rep(get_genedata("Prdm1_intron"), "TSS", "Prdm1_intron", limit = 0)
bf = bf[bf[!,:Timepoint].!= "30", :]
bf = bf[bf[!,:Timepoint].!= "90", :]
BF["Prdm1_intron"] = bf
bf = calculate_bf_by_rep(get_genedata("Il12b_intron"), "TSS", "Il12b_intron", limit = limitbs)
BF["Il12b_intron"] = bf
BFs = join_in_all_common_columns([BF[key] for key in keys(BF)]...)
CSV.write("BurstFrequencies_Genes.csv", BFs)
CSV.write("../SourceData/SupFig2_BF.csv", BFs)