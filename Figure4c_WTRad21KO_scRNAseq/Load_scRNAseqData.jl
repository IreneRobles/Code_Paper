function SingleCellExperiment.fit_single_cell_logistic_regression(singlecelldata; splitdataby = :Sample, assay = "CPM", fitparameter = "mu")
    cells = singlecelldata.colData[!,:CellID]
    b0 = []
    b1 = []
    p = Progress(length(cells))
    
    for cell in cells
        cell_index = find_cell_index(cell, singlecelldata)
        condition = find_cell_condition(cell, splitdataby, singlecelldata)
        genes__mus_in_condition = Array{Float64, 2}(singlecelldata.rowData[:, [Symbol(string(condition, "__", fitparameter))]])
         expression_cell_bool = Array{Float64,1}([i > 0 for i in singlecelldata.assays[assay][:, cell_index]])
        # fit each cell to a logistic regression
        model_ = LogisticRegression(penalty = "l2").fit(genes__mus_in_condition, expression_cell_bool)
        next!(p)
        push!(b0, model_.intercept_[1])
        push!(b1, model_.coef_[1])
    end
    
    singlecelldata.colData[!,:B0] = b0
    singlecelldata.colData[!,:B1] = b1
    
    return singlecelldata
    
end

# Perform the log(CPM + 1) transformation (Counts per million reads in each cell)
# Load the data
cd = SergiSingleCell_counts()
coldata = SergiSingleCell_colData();
coldata[!,:RowName] = [replace(ii, "."=> "_") for ii in coldata[!,:RowName]]

for ii in names(cd)
    stri = Symbol(replace("$ii", "."=> "_"))
    rename!(cd, ii => stri)
end
# Object that store Single cell Data
sce = SingleCellExp(cd, coldata)
sce = SingleCellExperiment.get_cells_with_this_characteristics(["WT", "RAD21"], :Genotype, sce)
sce = SingleCellExperiment.select_expressed_genes(sce, min_cells_expressing_gene = 20)
sce = SingleCellExperiment.cpm_transform(sce)
sce = SingleCellExperiment.ln_cpm_plus1_transform(sce)