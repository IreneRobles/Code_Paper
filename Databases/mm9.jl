module mm9

using DataFrames
using CSV
using NoLongerProblems
include("ProcessedData_mm9.jl")

function EnsemblIDtoGeneSymbol_table()
    ENSGeneID = string(ENV["Code"],"/Databases2/mm9/ENSGeneID_Symbol.txt")
    # This did not work in Julia 1.4
    #ENSGeneID = DataFrames.DataFrame!(CSV.File(ENSGeneID, delim = '\t', skipstart = 1, header = false))[2:end, [:Column1, :Column2]]
    
    ENSGeneID = DataFrames.DataFrame!(CSV.File(ENSGeneID, delim = '\t', skipto = 1, header = false))[2:end, [:Column1, :Column2]]
    rename!(ENSGeneID, :Column1 => :EnsemblID)
    rename!(ENSGeneID, :Column2 => :GeneSymbol)
end

function GeneLengths(; add_gene_symbols = false)
    genelen = string(ENV["Code"],"/Databases2/mm9/gene_length.txt")
    genelen = CSV.read(genelen, delim = '\t' ,skipto = 2, header = false)[:, [2, 3]]
    rename!(genelen, :Column2 => :EnsemblID)
    rename!(genelen, :Column3 => :GeneLength)
            
    if add_gene_symbols == true
        tb_symbols = EnsemblIDtoGeneSymbol_table()
        genelen = join(genelen, tb_symbols, on = :EnsemblID)
    
    end
            
    return genelen
      
end

export ensemblid_genetype_dict, ensemblid_genesymbol_dict, genesymbol_ensemblid_dict

function ensemblid_genetype_dict()
    tb = ProcessedData_mm9.mm9_genes_biomart()
    Dict(tb[:ensembl_gene_id], tb[:gene_biotype])
end


function ensemblid_genesymbol_dict()
    Dict(mm9.EnsemblIDtoGeneSymbol_table()[:EnsemblID],mm9.EnsemblIDtoGeneSymbol_table()[:GeneSymbol] )
end
        
function genesymbol_ensemblid_dict()
    Dict(mm9.EnsemblIDtoGeneSymbol_table()[:GeneSymbol],mm9.EnsemblIDtoGeneSymbol_table()[:EnsemblID] )
end
        


end