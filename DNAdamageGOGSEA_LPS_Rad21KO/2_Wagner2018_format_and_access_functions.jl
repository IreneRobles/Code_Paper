
function Wagner2018_DESEQ(inhibitor, dir = pwd())
    files = readdir(dir)
    filesbool = [occursin(".csv", ii) for ii in files]; files = files[filesbool]
    filesbool = [occursin(inhibitor, ii) for ii in files]; files = files[filesbool]
    filesbool = [startswith(ii,  "DESeq2_") for ii in files]; files = files[filesbool]
    if length(files) != 1 
        l = length(filesbool)
        ErrorException("$l files found")
    end
    
    table = dropmissing(CSV.read(files[1], DataFrames.DataFrame, missingstring = "NA"))
    table[!,:GeneSymbol] = translate_to_mouse(table[!, 1]; symbolstotranslate = "Gene name", translatedsymbols = "Mouse gene name")
    table[!,:EnsemblID] = translate_to_mouse(table[!, 1]; symbolstotranslate = "Gene name", translatedsymbols = "Mouse gene stable ID")
    
    return dropmissing(table)
end

function translate_to_mouse(genelist; symbolstotranslate = "Gene name", translatedsymbols = "Mouse gene name")
    ort = CSV.read(normpath(ENV["Code"], "Databases/GeneSets/mart_mouseortologstohuman.csv"), DataFrames.DataFrame)
    #Only one to one ortologs
    ort = ort[ort[!,"Mouse homology type"].=="ortholog_one2one", :]
    genedict = Dict(ort[!,symbolstotranslate], ort[!,translatedsymbols] )
    list =  [try genedict[split(ii, ".")[1]] catch; missing end for ii in genelist]
    return list
end

function entrezid(genelist; symbolstotranslate = "ensembl_gene_id", translatedsymbols = "entrezgene_id")
    entrez = CSV.read("../Code/Databases/GeneSets/Human_EnsemblID_.csv", DataFrames.DataFrame, delim = " ")
    genedict = Dict(ort[!,symbolstotranslate], ort[!,translatedsymbols] )
    list =  [try genedict[split(ii, ".")[1]] catch; missing end for ii in genelist]
    return list
end

function isbhatt(df)
    bhatt = Bhatt2012.inducible_genes_figure3()[!,:GeneSymbol]
    df[!,:IsBhatt] = [in(ii, bhatt) for ii in df[!,:GeneSymbol]]
    df
end

function updownnd(tb, suff1, suff2)


column1 = []
column2 = []

for ii in 1:nrow(tb)
    datum =  tb[ii, :]
    
    fc1 = "log2FoldChange_"*suff1
    padj1 = "padj_"*suff1
    d_fc1 =  datum[fc1]; d_padj1 =  datum[padj1]
        if d_padj1  < 0.05
            if d_fc1 > 0
                push!(column1, "up")
            else
                push!(column1, "down")
            end
        else
            push!(column1, "ND")
        end
    
     fc1 = "log2FoldChange_"*suff2
    padj1 = "padj_"*suff2
   
    d_fc1 =  datum[fc1]; d_padj1 =  datum[padj1]
        if d_padj1 < 0.05
            if d_fc1 > 0
                push!(column2, "up")
            else
                push!(column2, "down")
            end
        else
            push!(column2, "ND")
        end
end
    
tb[!,suff1] = column1

tb[!,suff2] = column2
tb[!,suff1*"_"*suff2] = suff1.*"_" .* column1 .*" ".*suff2.*"_" .*column2
tb[!,suff1*"_"*suff2]
tb
end

function printDEupDOWN(de; return_= false)
    n = nrow(de)
    dereg = de[de[!,:padj].<0.05, :]; n_dereg = nrow(dereg); fdereg = round(n_dereg/n*100, digits = 2)
    up = dereg[dereg[!,:log2FoldChange].>0, :];n_up = nrow(up); fup = round(n_up/n_dereg*100, digits = 2)
    down = dereg[dereg[!,:log2FoldChange].<0, :]; n_down = nrow(down);fdown = round(n_down/n_dereg*100, digits = 2)
    print("""
            DE genes = $n_dereg ($fdereg %)
                up genes = $n_up ($fup %)
                down genes = $n_down ($fdown %)
    """)
    
    if return_ == true
        return Dict("DE"=>dereg, "up"=>up, "down"=>down, "ND"=>de[de[!,:padj].>0.05, :])
    end
end


function geneOverlaps(listA, listB, GenomeSize)
a = R"""
library(GeneOverlap)
go.obj <- newGeneOverlap($listA,$listB,genome.size=$GenomeSize)
go.obj <- testGeneOverlap(go.obj)
print(go.obj)
a <-go.obj
"""
pval = round(rcopy(R"getPval(a)"), sigdigits = 2)
log2odds = round(log2(rcopy(R"getOddsRatio(a)")), sigdigits = 2)
return Dict("pvalue"=>pval,"log2OddsRatio"=>log2odds )
end

function venn2(a, b; kwargs...)
    
    n_all = length(intersect(a, b))

    n_1no = length(a) - n_all
    n_2no = length(b)  - n_all
    
    matplotlib_venn.venn2(subsets = (n_1no, n_2no, n_all); kwargs...)
end
        
function venn3(a, b, c; kwargs...)
    
    n_all = length(intersect(a, b, c))
    n_sim12 = length(intersect(a, b)) - n_all
    n_sim13 = length(intersect(a, c)) - n_all
    n_sim23 = length(intersect(b, c)) - n_all

    n_1no = length(a) - n_sim12 - n_sim13 - n_all
    n_2no = length(b) - n_sim12 - n_sim23 - n_all
    n_3no = length(c) - n_sim13 - n_sim23 - n_all

    matplotlib_venn.venn3(subsets = (n_1no, n_2no, n_sim12, n_3no, n_sim13, n_sim23, n_all); kwargs...)
end


function plotGeneOverlaps(A, B, n; title_ = "", labels = ["", ""])

venn2(A, B, set_labels = labels)
geneov = geneOverlaps(A, B, n)
pval = geneov["pvalue"]
log2odd = geneov["log2OddsRatio"]
title(title_)
annotate("""
    pvalue = $pval
    log2 Odds Ratio = $log2odd
    """, xy = (-0.7,0.5))
end

        

