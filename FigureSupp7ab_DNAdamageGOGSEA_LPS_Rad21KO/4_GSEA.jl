using DataFrames, CSV
using MSigDB
using RCall


function do_Hallmarks_GSEA(deseqtable; mouse_ensemblID = :ensembl_gene_id)
    
    deseqtable[!,:ensembl_gene_id] = deseqtable[!,mouse_ensemblID]
    deseqtable = dropmissing(deseqtable)
    codepath = ENV["Code"]
    
    
    
     R"""
    
library(devtools)
library(repr)
library(data.table)
library(fgsea)
library(ggplot2)
library("qusage")
library(dplyr)


tb1 <- read.csv(paste0($codepath, "/Databases/GeneSets/Human_EnsemblID_.csv"), sep = " ")
tb2 <- read.csv(paste0($codepath, "/Databases/GeneSets/mart_mouseortologstohuman.csv"))
tb2$ensembl_gene_id <- tb2$Gene.stable.ID

mg <- merge(tb1, tb2, by = "ensembl_gene_id")
mg$ensembl_gene_id <- mg$Mouse.gene.stable.ID
mg <- mg[mg$Mouse.orthology.confidence..0.low..1.high == 1, ]
mg <- mg[!duplicated(mg$ensembl_gene_id), ]

gseaDat <- $deseqtable
gseaDat <- merge(gseaDat, mg, by = "ensembl_gene_id")

gseaDat <- na.omit(gseaDat)
gseaDat <- gseaDat[!sapply(gseaDat$log2FoldChange, is.null), ]
ranks <- gseaDat$log2FoldChange
names(ranks) <- gseaDat$entrezgene_id

load(paste0($codepath, "/Databases3/MSigDB/human_H_v5p2.rdata"))


pathwaysH <- read.gmt(paste0($codepath, "/Databases3/MSigDB/h.all.v6.2.entrez.gmt"))
spli<- function(x){
    a = sub("HALLMARK_", "", x)
    return(a)
}


names(pathwaysH) <- lapply(names(pathwaysH), spli)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=50, maxSize = 10000, nperm=10000)
    """
    
    RCall.@rget fgseaRes
end

function joinGSEASforTable(gseas...; names = [string("x", ii) for ii in 1:length(gseas)])
    df = gseas[1][:, [:pathway, :NES]]
    rename!(df, :NES => Symbol(names[1]))
    for ii in 2:length(gseas)
        df2 = gseas[ii][:, [:pathway, :NES]]
        rename!(df2, :NES => Symbol(names[ii]))
        df = innerjoin(df, df2, on = :pathway)
    end
    df
end

function joinGSEASforTable_padj(gseas...; names = [string("x", ii) for ii in 1:length(gseas)])
    df = gseas[1][:, [:pathway, :padj]]
    rename!(df, :padj => Symbol(names[1]))
    for ii in 2:length(gseas)
        df2 = gseas[ii][:, [:pathway, :padj]]
        rename!(df2, :padj => Symbol(names[ii]))
        df = innerjoin(df, df2, on = :pathway)
    end
    df
end