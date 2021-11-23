
using DataFrames,CSV,Seaborn,Statistics
using NoLongerProblems, PrettyPlotting, NoLongerProblems_Pandas
using PyCall, RCall
import Pandas 

include(ENV["Code"]*"/../Code_Paper/Databases/Cuartero2018.jl")
include(ENV["Code"]*"/../Code_Paper/Databases/Bhatt2012.jl")

ind = Bhatt2012.inducible_genes_figure3()[!,[:GeneSymbol]]
ind[!,:Class] = ["Inducible" for ii in 1:nrow(ind)]

deseq = Cuartero2018.Cuartero2018Deseq("WT8_Minus_WTUT")
repress = deseq[deseq[!,:log2FoldChange] .< 0, :]
repress = repress[repress[!,:padj] .< 0.05, [:GeneSymbol]]
repress = unique(repress)
repress[!,:Class] = ["Repressed" for ii in 1:nrow(repress)]
repress

deseq = Cuartero2018.Cuartero2018Deseq("WT8_Minus_WTUT")

ctve = deseq[abs.(deseq[!,:log2FoldChange]) .< 0.5, :]
ctve1 = deseq[abs.(deseq[!,:padj]) .> 0.05, [:GeneSymbol]]
deseq = Cuartero2018.Cuartero2018Deseq("WT2_Minus_WTUT")

ctve = deseq[abs.(deseq[!,:log2FoldChange]) .< 0.5, :]
ctve2 = deseq[abs.(deseq[!,:padj]) .> 0.05, [:GeneSymbol]]

ctve = DataFrames.DataFrame(:GeneSymbol => intersect(ctve1[!,:GeneSymbol], ctve2[!,:GeneSymbol]))

ctve[!,:Class] = ["Constitutive" for ii in 1:nrow(ctve)]
if in("tables", readdir())==false
    mkdir("tables")
end
CSV.write("tables/MF_GeneClass.csv",join_in_all_common_columns(ind, repress, ctve))