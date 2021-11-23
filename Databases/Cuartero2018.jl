module Cuartero2018

using DataFrames
using CSV
using NoLongerProblems
using Statistics

include("../Databases/ProcessedData_mm9.jl")
include("../Databases/mm9.jl")

function get_FPKM()
    file = string("../Databases/Cuartero2018/FPKM_MF.csv")
     data = CSV.read(file, DataFrames.DataFrame)
    
    rename!(data, names(data)[1]=> :EnsemblID)
    
     data[!,:EnsemblID] =  Array{String, 1}(data[!,:EnsemblID])
    genenames = ProcessedData_mm9.one_gene_one_locus_databasetable()[!, [:GeneSymbol, :EnsemblID]]
     genenames[!,:EnsemblID] =  Array{String, 1}(genenames[!,:EnsemblID])    
    data = innerjoin(genenames, data, on=:EnsemblID)
end

function get_mean_FPKM()
    data_replicas = get_FPKM()
    cols = names(data_replicas)
    groups = unique([string(col)[5:end] for col in cols][3:end])

    groupdict = Dict()
    datadict = Dict()

    for group in groups
        groupdict[group] = Array{Float64, 1}()
        datadict[group] = select_columns_that_contain(data_replicas, group)
    end

    for i in 1:size(data_replicas)[1]
        for group in groups
            FPKM_mean =  Statistics.mean([jj for jj in datadict[group][i, :]])
            push!(groupdict[group], FPKM_mean)
        end
    end

    for group in groups
        data_replicas[!,Symbol(group)] = groupdict[group]
    end

    symbol_groups = [Symbol(group) for group in groups]
    push!(symbol_groups, :GeneSymbol)

    data_to_return = data_replicas[:, symbol_groups]

    rename!(data_to_return, :FLUT => :FL_0h)
    rename!(data_to_return, :FL2=> :FL_2h)
    rename!(data_to_return, :FL8=> :FL_8h)
    rename!(data_to_return, :WTUT=> :WT_0h)
    rename!(data_to_return, :WT2=> :WT_2h)
    rename!(data_to_return, :WT8=>:WT_8h)

    rename!(data_to_return, :FLi2=> :FL_i2h)
    rename!(data_to_return, :FLi8=> :FL_i8h)
    rename!(data_to_return, :WTi2=> :WT_i2h)
    rename!(data_to_return, :WTi8=> :WT_i8h)

    data_to_return = data_to_return[!, [:GeneSymbol, :WT_0h, :WT_2h, :WT_8h, :FL_0h, :FL_2h, :FL_8h, :WT_i2h, :WT_i8h, :FL_i2h, :FL_i8h]]
end




function Cuartero2018Deseq(name)
    yifandeseqfolder = ENV["Code"]*"/../Code_Paper/Databases/Cuartero2018/DE_result_norm2Rep_n_Spikes/"
    files = readdir(yifandeseqfolder)
     bool = [occursin(name, ii) for ii in files]
    file = files[bool]
     if length(file) == 1
        file = yifandeseqfolder*file[1]
        d = dropmissing(CSV.read(file, DataFrames.DataFrame, missingstring = "NA"))
        return rename!(d, :ensembl_id => :EnsemblID, :mgi_symbol => :GeneSymbol)
    else
        return files
    end
end

function GroseqDeseq(name)
    folder = ENV["Code"]*"/../Code_Paper/Databases/Cuartero2018/GROseq/"

    files = readdir(folder)
     bool = [occursin(name*"__DESeq2_AllGeneBody.csv", ii) for ii in files]
    file = files[bool]
     if length(file) == 1
        file = folder*file[1]
         d = innerjoin(dropmissing(CSV.read(file, DataFrames.DataFrame, missingstring = "NA")), mm9.EnsemblIDtoGeneSymbol_table(), on = :EnsemblID)
        return  d
    else
        return length(file)
    end
end

end