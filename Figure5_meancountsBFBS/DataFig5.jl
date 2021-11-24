ENV["Code"] = "../../Code"
for folder in readdir(ENV["Code"]); push!(LOAD_PATH, normpath(ENV["Code"], folder));end

include(ENV["Code"]*"/../Code_Paper/Code/meanmRNAcounts_BSBF.jl")

ifnb1 = CSV.read(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"CompleteSets/Ifnb1_mature.csv", DataFrames.DataFrame)
bf_ifnb1 = calculate_bf_by_rep(ifnb1, "TSS", "Ifnb1_mature")

il12b = CSV.read(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"CompleteSets/Il12b_mature_nascent.csv", DataFrames.DataFrame)
bf_il12b = calculate_bf_by_rep(il12b, "TSS", "Il12b")

peli1 = CSV.read(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"CompleteSets/Peli1_mature_nascent.csv", DataFrames.DataFrame)
bf_peli1 = calculate_bf_by_rep(peli1, "TSS", "Peli1")
                
Cxcl10 = CSV.read(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"CompleteSets/Cxcl10_mature.csv", DataFrames.DataFrame)
bf_Cxcl10 = calculate_bf_by_rep(Cxcl10, "TSS", "Cxcl10")

ifit1 = CSV.read(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"CompleteSets/Ifit1_mature.csv", DataFrames.DataFrame)
bf_ifit1 = calculate_bf_by_rep(ifit1, "TSS", "Ifit1")
   
                
bfs = join_in_all_common_columns([bf_ifnb1, bf_il12b, bf_ifit1, bf_peli1, bf_Cxcl10])
bfs[!,"log2 BF"] = log2.(bfs[!,"BF_TSS"])
bfs[!,"log2 BS"] = log2.(bfs[!,"BS_mean_TSS"])
bfs[!,"log2 mean mRNA counts"] = log2.(bfs[!,"MeanCountsTSS"]);

CSV.write("meanmRNAcounts_BSBF_v1.csv",bfs)


bfs = join_in_all_common_columns([bf_ifnb1, bf_il12b, bf_ifit1, bf_Cxcl10])
bfs[!,"log2 BF"] = log2.(bfs[!,"BF_TSS"])
bfs[!,"log2 BS"] = log2.(bfs[!,"BS_mean_TSS"])
bfs[!,"log2 mean mRNA counts"] = log2.(bfs[!,"MeanCountsTSS"]);

rename!(bfs, :MeanCountsTSS => :MeanCounts)

CSV.write("meanmRNAcounts_BSBF_v2.csv",bfs)