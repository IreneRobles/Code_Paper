function plot_gene_enh_correlation(tb_enh, tb_gene, genes; label = "WT", c = "blue", xy = (4, 0), enh = [])
    subtb_gene = tb_gene[!,[:EnsemblID, :log2FoldChange_Gene, :padj_Gene ]]

    jointb = innerjoin(tb_enh,subtb_gene, on = :EnsemblID)
    # Table with inducible genes and associated enhancers
    subtb = jointb[[in(ii, genes) for ii in jointb[!, :GeneSymbol_nearest]], :]
    if enh == []
    subtb_gene = subtb_gene[subtb_gene[!,:padj_Gene].< 0.05, :]
    jointb = innerjoin(tb_enh,subtb_gene, on = :EnsemblID)
    # Table with inducible genes and associated enhancers
    subtb = jointb[[in(ii, genes) for ii in jointb[!, :GeneSymbol_nearest]], :]
    # No Intragenic enhancers
    subtb = subtb[subtb[!,:Inter_Intragenic].=="Intergenic", :]
    # REmove nonDE enhancers
    subtb = subtb[subtb[!,:padj_Enh].<0.05, :]
    # Only inducible enhancers
    subtb = subtb[subtb[!,:log2FoldChange_Enh].>0, :]
    else
        subtb = subtb[[in(ii, enh) for ii in subtb[!, :coordinates]], :]
    end
    
    genes = unique(subtb[!,[:GeneSymbol_nearest, :log2FoldChange_Gene]])
    
    genes[!,:log2FoldChange_Enh] = [Statistics.mean(subtb[subtb[!,:GeneSymbol_nearest].==ii, :log2FoldChange_Enh]) for ii in genes[!,:GeneSymbol_nearest]]
    x = genes[!,:log2FoldChange_Enh]
    y = genes[!,:log2FoldChange_Gene]
    pval = round(R"""cor.test($x, $y)"""[3][1], sigdigits = 2)
    corr = round(R"""cor.test($x, $y)"""[4][1], sigdigits = 2)
    
    scatter(x, y, label = label, c = c)
    
    py"""
    import seaborn
    seaborn.regplot(x = $x, y = $y, scatter = 0, color = $c)
    """
    annotate("""P = $pval
        r = $corr""", xy = xy, c = c)

    ylabel("log2FoldChange Gene")
    xlabel("log2FoldChange Enhancers (avg)")
    pretty_axes2()
    squareplot()
    
    println(size(genes))
    
      return genes[!,:GeneSymbol_nearest], subtb[!,:coordinates], genes
end

function maketablesforanalysis(gro_gene, gro_enh, filename, genotype; contacts = false)
gro_gene =  unique(gro_gene)
mm9table = ProcessedData_mm9.one_gene_one_locus_databasetable()[!,[:EnsemblID]]
gro_gene = innerjoin(mm9table,gro_gene, on = :EnsemblID)
mm9table = ProcessedData_mm9.one_gene_one_locus_databasetable()[!,[:EnsemblID, :GeneSymbol, :chr, :start, :end, :tss]]
mm9table = innerjoin(mm9table, DataFrames.DataFrame(GeneSymbol = geneset), on = :GeneSymbol)
tads = ""
    if !contacts
    tads = bed_readerTADs(normpath(ENV["Code"]*"/../Code_Paper/Databases/GSE115524/", "TADs_list_25kb_iced_GSE115524_HiC_Macrophage.bed"))
    else
    tads = bed_readerTADs(normpath(ENV["Code"]*"/../Code_Paper/Databases/GSE115524/","Contact_Domains_GSE115524_HiC_Macrophage_R1R2_keepSmall.bed"))
    end

#Get the TADs per gene based on their TSS
TADgenes = []
for ii in 1:nrow(mm9table)
    chr = mm9table[ii, :chr]
    chr_tad = tads[tads[!,:TADchr].==string("chr",chr), :]
    tss = mm9table[ii, :tss]
    symbol = mm9table[ii, :GeneSymbol]
     gene_gro = gro_gene[gro_gene[!,:GeneSymbol].== symbol, :]
 
  try
    datum =  chr_tad[(chr_tad[!,:TADstart].< tss).* (tss .< chr_tad[!,:TADend]), 1:4]
    datum[!,:GeneSymbol] = [symbol]
     datum[!,:Genelog2FoldChange] = [mean(gene_gro[:, :log2FoldChange_Gene])]
    datum[!,:Genepadj] = [mean(gene_gro[:, :padj_Gene])]
    push!(TADgenes,datum)
    catch
 
    end
   
end

TADgenes  = join_in_all_common_columns(TADgenes...)

#Get the TADs per enhancer based on their middlepoint
TADenh = []
for ii in 1:nrow(gro_enh)
    chr = gro_enh[ii, :chr]
    chr_tad = tads[tads[!,:TADchr].==string("chr",chr), :]
   
    tss = mean([gro_enh[ii, :start], gro_enh[ii, :end]])
   try
        datum =  chr_tad[(chr_tad[!,:TADstart].< tss).* (tss .< chr_tad[!,:TADend]), [4]]
        datum[!,:EnhCoordinates] = [gro_enh[ii, :coordinates]]
        datum[!,:Enhlog2FoldChange] = [gro_enh[ii, :log2FoldChange_Enh]]
        datum[!,:Enhpadj] = [gro_enh[ii, :padj_Enh]]
    push!(TADenh,datum)
    catch
    end
end

TADenh = join_in_all_common_columns(TADenh...)

tb = innerjoin(TADgenes, TADenh, on = "TAD" )
tb[!,:Genotype] = [genotype for ii in 1:nrow(tb)]
CSV.write("$filename", tb)
    return tb
end

