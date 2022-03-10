using StatsBase
function show_distances2(genefolder,nascent,ehn,suff; limit = 4, limit_gene = 1, limit_enh =1) 
    tb = linked_data(genefolder,nascent,ehn,suff) 
    tb = tb[tb[!,:Timepoint].>0, :]
    if nascent == "Prdm1_intron"
        tb = tb[tb[!,:Timepoint].==60, :]
    end
    if nascent == "Il12b_intron" || nascent == "Egr2_intron"
        tb2 = linked_data(genefolder, nascent*"_BD",ehn*"_BD",suff) 
        tb2 = tb2[tb2[!,:Timepoint] .== "DMSO-LPS", :]
        tb = join_in_all_common_columns(tb, tb2)
        if nascent == "Il12b_intron"
            tb2 = linked_data(genefolder, nascent*"_timecourse",ehn*"_timecourse",suff) 
            tb2 = tb2[tb2[!,:Timepoint] .== 90, :]
            tb2[!,:Rep]= tb2[!,:Rep].+20
            tb = join_in_all_common_columns(tb, tb2)
        end
    end
    pairs = get_pairs_with_distances(tb)
    
    pairs = pairs[parse.(Float64,pairs[!,:locus1_Gene_Enh]).< limit, :]
    sort!(pairs, [:Genotype], rev = true)
    
    CSV.write("../SourceData/Fig7"*string("_",nascent, "_",ehn)*".csv", pairs)


    # generate the samples
    
     g1 = parse.(Float64,pairs[pairs[!,:Genotype].=="WT", :locus1_Gene_Enh])
    gcdf = ecdf(g1)
     
    g2 = parse.(Float64,pairs[pairs[!,:Genotype].=="Rad21KO", :locus1_Gene_Enh])
    gcdf2 = ecdf(g2)
    
    range = 0:0.005:limit
    plt.plot(range,[gcdf(x) for x in range],  label = "WT", c = "darkgray")
    plt.plot(range,[gcdf2(x) for x in range], label = "Rad21KO", c = "red")

    
    pretty_axes2()
   ylabel("Cumulative Probability")
    xlabel("Distance Gene-Enhancer (um)")
    line075black()

    squareplot()
    pval =  round(pvalue(HypothesisTests.ApproximateTwoSampleKSTest(g1,g2)), sigdigits = 2)

     title("$nascent - $ehn \nP=$pval")
   
    
end

function show_distances_BD(genefolder,nascent,ehn,suff; limit = 4, limit_gene = 1, limit_enh =1, subplots = [(1,2,1), (1,2,2)]) 
    tb = linked_data(genefolder,nascent,ehn,suff) 
    if genefolder == "Prdm1"
        tb = tb[tb[!,:Timepoint].==60, :]
    end
    pairs = get_pairs_with_distances(tb)
    
    pairs = pairs[parse.(Float64,pairs[!,:locus1_Gene_Enh]).< limit, :]
    sort!(pairs, [:Genotype, :Timepoint], rev = true)
    
    wt = pairs[pairs[!,:Genotype].=="WT", :]
    rad21 = pairs[pairs[!,:Genotype].=="Rad21KO", :]

    # generate the samples for WT
    
    subplot(subplots[1]...)
    
    g_dmso = parse.(Float64,wt[wt[!,:Timepoint].=="DMSO-LPS", :locus1_Gene_Enh])
    gcdf_dmso = ecdf(g_dmso)
    g_bd1 = parse.(Float64,wt[wt[!,:Timepoint].=="BD1-LPS", :locus1_Gene_Enh])
    gcdf_bd1 = ecdf(g_bd1)
    g_bd2 = parse.(Float64,wt[wt[!,:Timepoint].=="BD2-LPS", :locus1_Gene_Enh])
    gcdf_bd2 = ecdf(g_bd2)
    
    pval1 =  round(pvalue(HypothesisTests.ApproximateTwoSampleKSTest(g_dmso,g_bd1)), sigdigits = 2)
    pval2 =  round(pvalue(HypothesisTests.ApproximateTwoSampleKSTest(g_dmso,g_bd2)), sigdigits = 2)

    
    range = 0:0.005:limit
    plt.plot(range,[gcdf_dmso(x) for x in range],  label = "WT_DMSO-LPS", c = "black")
    plt.plot(range,[gcdf_bd1(x) for x in range], label = "WT_BD1-LPS P=$pval1", c = "#2b79d9")
    plt.plot(range,[gcdf_bd2(x) for x in range], label = "WT_BD2-LPS P=$pval2", c = "#6936e0")
    
    legend()
        pretty_axes2()
   ylabel("Cumulative Probability")
    xlabel("Distance Gene-Enhancer (um)")
    line075black()
    squareplot()
    title(nascent)

    subplot(subplots[2]...)
    g_dmso = parse.(Float64,rad21[rad21[!,:Timepoint].=="DMSO-LPS", :locus1_Gene_Enh])
    gcdf_dmso = ecdf(g_dmso)
    g_bd1 = parse.(Float64,rad21[rad21[!,:Timepoint].=="BD1-LPS", :locus1_Gene_Enh])
    gcdf_bd1 = ecdf(g_bd1)
    g_bd2 = parse.(Float64,rad21[rad21[!,:Timepoint].=="BD2-LPS", :locus1_Gene_Enh])
    gcdf_bd2 = ecdf(g_bd2)
    
    pval1 =  round(pvalue(HypothesisTests.ApproximateTwoSampleKSTest(g_dmso,g_bd1)), sigdigits = 2)
    pval2 =  round(pvalue(HypothesisTests.ApproximateTwoSampleKSTest(g_dmso,g_bd2)), sigdigits = 2)

    
    range = 0:0.005:limit
    plt.plot(range,[gcdf_dmso(x) for x in range],  label = "Rad21KO_DMSO-LPS", c = "red")
    plt.plot(range,[gcdf_bd1(x) for x in range], label = "Rad21KO_BD1-LPS P=$pval1", c = "#2b79d9")
    plt.plot(range,[gcdf_bd2(x) for x in range], label = "Rad21KO_BD2-LPS P=$pval2", c = "#6936e0")
     

    legend()
    title(nascent)
    pretty_axes2()
   ylabel("Cumulative Probability")
    xlabel("Distance Gene-Enhancer (um)")
    line075black()
    squareplot()

     #title("$nascent - $ehn \nP=$pval")
   plt.tight_layout()
    
end



function pair_summary(folder, nascent, enh, suff ; limit = 10)
    tb = linked_data(folder, nascent, enh, suff)
    
    tb = tb[tb[!,:Timepoint].>0, :]
    if nascent == "Prdm1_intron"
        tb = tb[tb[!,:Timepoint].==60, :]
    end
    if nascent == "Il12b_intron" || nascent == "Egr2_intron"
        tb2 = linked_data(folder, nascent*"_BD",enh*"_BD",suff) 
        tb2 = tb2[tb2[!,:Timepoint] .== "DMSO-LPS", :]
        tb2[!,:Rep]= tb2[!,:Rep].+5
        tb = join_in_all_common_columns(tb, tb2)
        if nascent == "Il12b_intron"
            tb2 = linked_data(folder, nascent*"_timecourse",enh*"_timecourse",suff) 
            tb2 = tb2[tb2[!,:Timepoint] .== 90, :]
            tb2[!,:Rep]= tb2[!,:Rep].+20
            tb = join_in_all_common_columns(tb, tb2)
        end
    end
     pairs = get_pairs_with_distances(tb)
    pairs = pairs[parse.(Float64,pairs[!,:locus1_Gene_Enh]).< limit, :]
    wt = pairs[pairs[!,:Genotype].=="WT", :]
    rad = pairs[pairs[!,:Genotype].=="Rad21KO", :]
    wttb = tb[tb[!,:Genotype].=="WT", :]
    radtb = tb[tb[!,:Genotype].=="Rad21KO", :]
    
    
    DataFrames.DataFrame(
        Intron = [nascent],
        Enhancer = [enh],
        Reps = [length(unique(tb[!,:Rep]))],
        N_Cells = [nrow(tb)],
        N_Cells_WT = [nrow(wttb)],
        N_Cells_Rad21KO = [nrow(radtb)],
        N_Pairs = [nrow(pairs)],
        N_Pairs_WT = [nrow(wt)],
        N_Pairs_Rad21KO = [nrow(rad)],
    )
end