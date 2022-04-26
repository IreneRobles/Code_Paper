ifnb1 = CSV.read(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"CompleteSets/Ifnb1_mature.csv", DataFrames.DataFrame)
bf_ifnb1 = bootstrap_mean_std(ifnb1, "Ifnb1")
bf_ifnb1 = bf_ifnb1[parse.(Int, bf_ifnb1[!,:Timepoint]).==90, :]


il12b = CSV.read(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"CompleteSets/Il12b_mature_nascent.csv", DataFrames.DataFrame)
bf_il12b = bootstrap_mean_std(il12b, "Il12b")
bf_il12b = bf_il12b[parse.(Int, bf_il12b[!,:Timepoint]).==90, :]
                
Cxcl10 = CSV.read(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"CompleteSets/Cxcl10_mature.csv", DataFrames.DataFrame)
bf_Cxcl10 = bootstrap_mean_std(Cxcl10, "Cxcl10")
bf_Cxcl10 = bf_Cxcl10[bf_Cxcl10[!,:Timepoint].!="0", :]


ifit1 = CSV.read(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"CompleteSets/Ifit1_mature.csv", DataFrames.DataFrame)
bf_ifit1 = bootstrap_mean_std(ifit1, "Ifit1")
bf_ifit1 = bf_ifit1[bf_ifit1[!,:Timepoint].=="180", :]

peli1 = CSV.read(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"CompleteSets/Peli1_mature_nascent.csv", DataFrames.DataFrame)
bf_peli1 = bootstrap_mean_std(peli1, "Peli1")
bf_peli1 = bf_peli1[parse.(Int, bf_peli1[!,:Timepoint]).==90, :]


bfs = join_in_all_common_columns([bf_peli1, bf_il12b,bf_ifit1,bf_Cxcl10,bf_ifnb1])
bfs = bfs[bfs[!,:Genotype].=="WT", :]