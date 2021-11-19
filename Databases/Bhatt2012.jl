module Bhatt2012

using DataFrames
using PrettyPlotting
using CSV
using PyPlot


export all_genes_figure2, inducible_genes_figure3, show_me_fractions

function all_genes_figure2()
    file = string("../Databases/Bhatt2012/Figure2.csv")
    data = readtable(file)
    rename!(data, :_Official_Gene_Name, :GeneSymbol)
    rename!(data, :NAME, :Gene_Refseq)
    
end

function isbhatt(df)
    bhatt = Bhatt2012.inducible_genes_figure3()[!,:GeneSymbol]
    df[!,:IsBhatt] = [in(ii, bhatt) for ii in df[!,:GeneSymbol]]
    df
end

function inducible_genes_figure3()
    file = string("../Databases/Bhatt2012/Figure3.csv")
    data = CSV.read(file, DataFrames.DataFrame)
    rename!(data, :Gene => :GeneSymbol)
end


#Old functions

"""
using RNAseq_v2, RNAseq_Sergi
function show_me_fractions(df::DataFrame, genes...)
    n = length(genes)
    count = 0
    for gene in genes
        if in(gene, df[!,:GeneSymbol])
            count+=1
            if length(genes) > 1
                subplot(n, 1, count)
            end
            datum = show_me_gene(df, gene)
            labels = [0, 0.25, 0.5, 1, 2]
            chromatin = [:chrom0, :chrom0_25, :chrom0_5, :chrom1, :chrom2]
            nucleoplasm = [:nuc0, :nuc0_25, :nuc0_5, :nuc1, :nuc2]
            cytoplasm = [:cyto0, :cyto0_25, :cyto0_5, :cyto1, :cyto2]
            values =  reshape(Array(datum[1,chromatin]),5)
            plot(labels, values, label = "chromatin", lw = 2)
            values =  reshape(Array(datum[1,nucleoplasm]),5)
            plot(labels, values, label = "nucleoplasm")
            values =  reshape(Array(datum[1,cytoplasm]),5)
            plot(labels, values, label = "cytoplasm")
            legend_out_of_plot()
            title(gene)
            pretty_axes()
            xticks(labels)
            ylabel("RPKM")
        end
    end
    xlabel("time(h after LPS treatment)")
    plt[:tight_layout]()
end
"""

end
