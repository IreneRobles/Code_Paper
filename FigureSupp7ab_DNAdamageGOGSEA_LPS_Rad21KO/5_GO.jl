rscript = R"""

library("goseq")

perform_goseq <- function(deseq, 
                          analysisname = "Test",
                          genome = "mm9", 
                          geneid = 'ensGene',
                          plot.fit = F,
                          test.cats=c("GO:BP", "GO:MF"),
                          upregulated = T # If upregulated == true, then the test is on upregulated genes, otherwise the test in on downregulated genes
                         ){
    
    degenes<-deseq$Test
    names(degenes) <- deseq$EnsemblID

    # remove duplicate gene names
    degenes<-degenes[match(unique(names(degenes)), names(degenes))]
   
    # We first need to obtain a weighting for each gene, depending on its length, given by the PWF
   
    pwf=nullp(degenes,genome=genome, geneid, plot.fit=plot.fit)
    

    
    go<-goseq(pwf,genome=genome,geneid, test.cats=test.cats)
    enriched.GO= p.adjust(go$over_represented_pvalue, method="BH")
    go$over_represented_padj <- enriched.GO
    
    write.csv(go, file = paste0(analysisname, "__goseq.csv"))
    
    return(paste0(analysisname,"__goseq.csv"))
    
}

"""



function perform_goseq(
        genelist;
        analysisname = "Analysis",
        genome = "mm9", 
        geneid = "ensGene", 
        plot_fit = false, 
        test_cats = ["GO:BP", "GO:MF","KEGG"], 
        upregulated = true)
    

     d = .! nonunique(genelist[:, [:EnsemblID]])
     d = genelist[d, :]
    
    R"""
    d = $d
    go <- perform_goseq(d, 
    upregulated = $upregulated, 
    analysisname = $analysisname, 
    test.cats = $test_cats)     
    """
    
    
end

 
function barplot_go(up_down; n = 10)
    
    file = up_down
    file[!,Symbol("-log10(overrepresented_padj)")] = [-log10(i) for i in file[!,:over_represented_padj]]
    
    pd = Pandas.DataFrame(file[1:n, :])
    
    barplot(x = "-log10(overrepresented_padj)", y="term", data=pd, palette="Blues_d")
    
    plot([-log10(0.05), -log10(0.05)], [ -0.5 , n - 0.5 ], linewidth=2, c = "red")
    
    annotate("padj = 0.05" , xy = [-log10(0.05),  n ], color = "darkred")
    
    ylim(-0.75 , n+1)
    
    axes = gca()
    axes[:spines]["top"][:set_visible](false) # Hide the top edge of the axis
    axes[:spines]["right"][:set_visible](false) # Hide the right edge of the axis
    axes[:xaxis][:set_ticks_position]("bottom")
    axes[:yaxis][:set_ticks_position]("left")
end



function plot_category(goterm)
    minitb = maxitb[maxitb[!,:category].==goterm, :]
    title_= minitb[1,:term]
    file = sort(minitb, :Test)
    n = 24
    file[!,Symbol("-log10(overrepresented_padj)")] = [-log10(i) for i in file[!,:over_represented_padj]]
    pd = Pandas.DataFrame(file[1:n, :])
    CSV.write("../SourceData/SupFig7_"*replace(title_, " "=>"_")*".csv",file[1:n, :])
    barplot(x = "-log10(overrepresented_padj)", y="Test", data=pd, palette="Blues_d")
    
    plot([-log10(0.05), -log10(0.05)], [ -0.5 , n - 0.5 ], linewidth=2, c = "red")
    
    annotate("padj = 0.05" , xy = [-log10(0.05),  n ], color = "darkred")
    
    ylim(-0.75 , n+1)
    
    PrettyPlotting.pretty_axes2()
    plt.title(title_)
    minitb
end

