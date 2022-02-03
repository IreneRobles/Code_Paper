col_colors_dict = Dict(
    "WT_UT" => "#f0f0f0",
    "WT_2H" => "#bdbdbd",
    "WT_8H" => "#636363",
    "RAD21_UT" => "#fee0d2",
    "RAD21_2H" => "#fc9272",
    "RAD21_8H" => "#de2d26",
    "RAD21" => "red",
    "WT" => "darkgrey",
    0 => "#deebf7",
    2 => "#9ecae1",
    8 => "#3182bd",
    ) 


sce.colData[!,"LPS"] = [replace(replace(replace(split(ii, "_")[2], "UT" => 0), "2H" => 2), "8H" => 8) for ii in sce.colData[!,"Sample"]]


subsce = select_these_genes(bhattgenes[!,:GeneID], sce)

# Bool for 2H LPS
bool1 = subsce.rowData[!,:WT_2H__alpha].>=alpha
# Bool for 8H LPS
bool2 = subsce.rowData[!,:WT_8H__alpha].>=alpha
# Bool for either in 2H or 8H
bool = (bool1.+bool2).>0
println(string("Percent ", alpha*100, "%"))
println(string("Total genes considered ", sum(bool)))

subsubsce = filter_genes(bool,subsce)
subsubsce = SingleCellExperiment.Shalek2014_module_score(collect(bhattgenes[!,:GeneID]), subsubsce,fitparameter = "mu", modulescore_name = :BhattGenesScore, untreated_pattern = "UT",comparedtothissample = "WT", assay = "CPM")
subsubsce = sort_cells!(subsubsce, cols = [ :LPS, :Genotype,:BhattGenesScore], rev = [false, true, false])



g = Seaborn.clustermap(subsubsce.assays["lnCPMplus1"], figsize = (7, 4),yticklabels=false,xticklabels=false, col_cluster = false, col_colors = [col_colors_dict[ii] for ii in subsubsce.colData[!,:Genotype]], cmap="coolwarm")

ax = g.ax_heatmap
ax.set_ylabel("LPS inducible genes")
ax.set_xlabel("Cells")
ax.tick_params(axis="both", which="both", length=0)

savefigwithtext("scRNAseq_heatmap_bhattgenes_WTRad21KO"*"percent"*string(alpha*100)*".pdf")