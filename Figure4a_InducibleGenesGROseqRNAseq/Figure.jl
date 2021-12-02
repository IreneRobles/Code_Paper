function figure_de_gro_cor_induciblegenes(de, gro, pref1, pref2)
    figure(figsize = (3, 3))

    tb = innerjoin(de, gro, on = :EnsemblID, makeunique = true)
    tb = innerjoin(dropmissing!(GeneClassTB),tb , on = :GeneSymbol)
    tb = tb[tb[!,:Class].=="Inducible", :]
    
    
    tbde1 = []
    tbde2 = []
    for ii in 1:nrow(tb)
        if tb[ii, "padj"] .> 0.05
            push!(tbde1, split(pref1)[1]*":ns")
        elseif tb[ii, "log2FoldChange"] .> 0
            push!(tbde1, split(pref1)[1]*":up")
        elseif tb[ii, "log2FoldChange"] .< 0
            push!(tbde1, split(pref1)[1]*":down")
        end
        if tb[ii, "padj_1"] .> 0.05
            push!(tbde2, split(pref2)[1]*":ns")
        elseif tb[ii, "log2FoldChange_1"] .> 0
            push!(tbde2, split(pref2)[1]*":up")
        elseif tb[ii, "log2FoldChange_1"] .< 0
            push!(tbde2, split(pref2)[1]*":down")
        end
    end
    
    tb[!,"DE1"] = tbde1
    tb[!,"DE2"] = tbde2
    tb[!,"DE"] = tbde1.*" ".*tbde2

    name  = replace(pref1, " "=>"")*"_"*replace(pref2, " "=>"")
    
    tb = tb[tb[!,"DE"].!=split(pref1)[1]*":ns".*" ".*split(pref2)[1]*":ns", :]
    tb = tb[tb[!,"DE1"].!=split(pref1)[1]*":ns", :]
    sort!(tb, "DE", rev = true)
    
    pd = Pandas.DataFrame(tb)
    
    Seaborn.lmplot(data = pd, x = "log2FoldChange", y = "log2FoldChange_1", hue = "DE", fit_reg = false, scatter_kws=Dict("alpha"=>0.5),
        #Matthias colours:
        #RNA-seq up GRO-seq up red, 
        #RNA-seq up GRO-seq NS light grey, RNA-seq up GRO-seq down orange
        #RNA-seq down GRO-seq down blue, 
        #RNA-seq down GRO-seq NS medium grey, 
        #RNA-seq down GRO-seq up turquoise 
        
        
        palette = Dict(
        "RNAseq:up GROseq:up"=> "red",
        "RNAseq:up GROseq:ns"=> "lightgrey",
        "RNAseq:up GROseq:down"=> "orange",
        "RNAseq:down GROseq:down"=> "blue",
        "RNAseq:down GROseq:ns"=> "darkgrey",
        "RNAseq:down GROseq:up"=> "turquoise",
        
        )
    
    )
    
    xlabel("log2FC \n"*pref1)
    ylabel("log2FC \n"*pref2)
    legend_removal()
    x = tb[!,"log2FoldChange"]; y = tb[!,"log2FoldChange_1"]
    print(R"cor.test($x,$y)")
    squareplot()
    fname = name*"_LPSinducible.svg"
    savefigwithtext(fname)
end



