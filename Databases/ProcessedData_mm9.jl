module ProcessedData_mm9

using DataFrames
using CSV

export mm9_genes_biomart, one_gene_one_locus_databasetable, mm9_genes_biomart_standard_name

function mm9_genes_biomart()
    file = string("../Databases/ProcessedData_mm9/mm9_biomart_genes.csv")
    database = CSV.read(file, DataFrames.DataFrame)
    database
end

function mm9_genes_biomart_standard_name()
    file = string("../Databases/ProcessedData_mm9/mm9_biomart_genes.csv")
    database = readtable(file)
    rename!(database, :ensembl_gene_id => :EnsemblID)
    rename!(database, :external_gene_id => :GeneSymbol)
    rename!(database, :chromosome_name => :chr)
    rename!(database, :gene_biotype => :Gene_type)
    rename!(database, :end_position => :end)
    rename!(database, :start_position => :start)
    database
end

function one_gene_one_locus_databasetable(; kind = "longest_locus")
    if kind == "longest_locus"
        file = string("..", "/Databases/ProcessedData_mm9/genedatabasefornetworkskeleton_kindlongestlocuspergene_mm9.csv")
        database = CSV.read(file, DataFrames.DataFrame)
        try rename!(database, :_end => :end) catch end
        return database
    else
        println("This kind is not available")
    end
end
end
