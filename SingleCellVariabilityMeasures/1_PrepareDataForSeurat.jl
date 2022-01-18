ENV["Code"] = "../../Code"
for folder in readdir(ENV["Code"]); push!(LOAD_PATH, normpath(ENV["Code"], folder));end

using SergiSingleCell2
using SingleCellExperiment
using DataFrames
using CSV

cd = SergiSingleCell2.SergiSingleCell_counts()
coldata = SergiSingleCell2.SergiSingleCell_colData();
coldata[!,:RowName] = [replace(ii, "."=> "") for ii in coldata[!,:RowName]]
coldata[!,:RowName] = [replace(ii, "_"=> "") for ii in coldata[!,:RowName]]
coldata[!,:RowNameName] = coldata[!,:RowName].*"_".*coldata[!,:Genotype].*coldata[!,:Timepoint]

for ii in names(cd)
    stri = Symbol(replace("$ii", "."=> ""))
    stri = Symbol(replace("$ii", "_"=> ""))
    global cd = rename!(cd, ii => stri)
end

sce = SingleCellExp(cd, coldata)
sce = SingleCellExperiment.get_cells_with_this_characteristics(["WT", "RAD21"], :Genotype, sce)

function prepare_data_seurat(scexp; dir = "SergiSeurat")
    if !in(dir, readdir()) mkdir(dir) end 
    
    CSV.write(dir*"/counts.csv",  DataFrame(scexp.counts), header=false)
    CSV.write(dir*"/metadata.csv",  DataFrame(scexp.colData))
    CSV.write(dir*"/genes.csv",  DataFrame(A = scexp.rowData[!,:GeneID]), header=false)
    CSV.write(dir*"/cells.csv",  DataFrame(A = scexp.colData[!,:CellID]), header=false)
end

prepare_data_seurat(sce)