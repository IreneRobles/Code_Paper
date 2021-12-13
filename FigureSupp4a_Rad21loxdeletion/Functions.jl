function NoLongerProblems.rename_column_and_merge_dfs(symbol, dict)
    Keys = collect(keys(dict))
    values = [dict[i] for i in Keys]
    
    NoLongerProblems.rename_column_and_merge_dfs(symbol, values, Keys)
    
end

function get_mean_Cq_by_target(df::DataFrame)
    df[!,:Sample_Target] = [string(df[i, :Sample], "__", df[i, :Target]) for i in 1:nrow(df)]
     sample_targets = unique(df[!,:Sample_Target])
    
    new_df = DataFrame()
    new_df[!,:Sample_Target] = sample_targets
    new_df[!,:Sample] = [string(split(i, "__")[1]) for i in sample_targets]
    new_df[!,:Target] = [split(i, "__")[end] for i in sample_targets]
   # new_df[:Genotype] = [split(i, "_")[1] for i in sample_targets]
    new_df[!,:Mean_Cq] = [mean(df[df[!,:Sample_Target] .== s, :Cq])  for s in sample_targets]
    new_df[!,:Std_Cq] = [StatsBase.std(df[df[!,:Sample_Target] .== s, :Cq])  for s in sample_targets]
    new_df
end

function get_delta_Cq_by_target(df::DataFrame; housekeeping = ["Actb", "Hprt"])
     mean_cq_df = get_mean_Cq_by_target(df)
    samples = unique(mean_cq_df[!,:Sample])
    targets = unique(mean_cq_df[!,:Target])
    
    new_df = DataFrame()
    new_df[!,:Samples] = samples

    averages_cq = []
    
    for sam in samples
        # Get the mean Cq value for the housekeeping genes in each sample
        data_hk =  mean_cq_df[.&((mean_cq_df[!,:Sample] .== sam), [in(i, housekeeping) for i in mean_cq_df[!,:Target]]), :]
        push!(averages_cq, StatsBase.mean(data_hk[!,:Mean_Cq]))
    end
    
     new_df[!,:Norm_Cq] = averages_cq
    
    
    
    for target in targets
        print(target)
        Cq_col_name = Symbol(string("Cq_", target))
        dCq_col_name = Symbol(string("dCq_", target))
        
        Cq_values = []
        dCq_values = []
        
        for sam in samples
                
            Cq = 0
            dCq = 0
                
            try
        
            Cq = mean_cq_df[.&(
                    (mean_cq_df[!,:Sample] .== sam),
                    (mean_cq_df[!,:Target] .== target)                    
                    ), :Mean_Cq][1]
             
            dCq = new_df[new_df[!,:Samples] .== sam, :Norm_Cq][1]
                catch;
                    
                end
            
                
            push!(Cq_values, Cq)
            push!(dCq_values, Cq - dCq)
                
            
        end
         Cq_values
        new_df[!,Cq_col_name] = Cq_values
        new_df[!,dCq_col_name] = dCq_values
    end
    
    new_df
    
    
    
end




function get_deltadelta_Cq(df_withdCq, control_sample)
    
    col_n = names(df_withdCq)
    dCqcolumns = NoLongerProblems.columns_containing(df_withdCq, "dCq")
    
    columns = append!(["Samples"], col_n)
        
    new_df = df_withdCq[!, unique(columns)]
    
    samples = new_df[!,:Samples]
    
    control_condition = new_df[new_df[!,:Samples] .== control_sample, :]
    
    
    for dCq in dCqcolumns
        
        ddCq = Symbol(string("d", dCq))
        ddCq_values = []
        expfolf_values = []
        
        for sam in samples
        
           
            dcq = df_withdCq[df_withdCq[!,:Samples] .== sam, dCq][1, 1]
             
            dcq_control = control_condition[1, dCq]
            
            ddcq = dcq - dcq_control
            
            
            push!(ddCq_values, ddcq)
            push!(expfolf_values, 2^(-ddcq))
        end
        
        new_df[!,ddCq] = ddCq_values
        new_df[!,Symbol(string(string(dCq)[5:end], "_FoldChange"))] = expfolf_values
    
        
    end
        try
    new_df[!,:Genotype] = [split(i, "_")[1] for i in samples]
    new_df[!,:Timepoint] = [split(i, "_")[2] for i in samples]
            catch; 
        end
    new_df
    
    
end
    

function final_figure_rad21(ddf; 
        col = :Ctcf_NormLevel, 
        control = "WT", 
        nametable = "rep2rep3.csv",
        cs = ["blue", "orange", "green", "red"]
    )
     

 

    ddf[!,:Genotype] =[if ii[end] == '5' "Rad21KO_5Days" else ii end for ii in ddf[!,:Genotype] ] 
    ddf[!,:Genotype] =[if ii[end] == '6' "Rad21KO_6Days" else ii end for ii in ddf[!,:Genotype] ] 
    ddf[!,:Genotype] =[if ii[end] == 'O' "Rad21KO_7Days" else ii end for ii in ddf[!,:Genotype] ] 
    
     ddf = ddf[ddf[!,:Day].> 4, :]

    
    sort!(ddf, :Genotype, rev = true)
    
    gens = ["Rad21lox", "Rad21KO_5Days",  "Rad21KO_6Days", "Rad21KO_7Days"]

    

     pd = Pandas.DataFrame(ddf)
    
    pd[string(col)] = [100 - ii*100 for ii in pd[string(col)]]
    
    
     Seaborn.barplot(x="Genotype", y=string(col), data=pd, palette =  "Reds", order = ["Rad21lox", "Rad21KO_5Days","Rad21KO_6Days", "Rad21KO_7Days"],)
    
    xs = collect(0:length(unique(ddf[!,:Genotype]))-1)
    ys = [maximum(pd[string(col)])+maximum(pd[string(col)])*0.2*i for i in xs]
    
    for i in 2:length(xs)
        plot([0, xs[i]], [ys[i], ys[i]], linewidth=1, color="black")
        # Asses statistical significance
        f(x, gen) = [i == gen for i in x[!,:Genotype]]
        wt = [Float64(i) for i in collect(ddf[f(ddf, control), col])]
        ko = [Float64(i) for i in collect(ddf[f(ddf, gens[i]), col])]
        test = scipy.stats.f_oneway(wt, ko)
        pval = test[2]

        if pval <= 0.0001
            tag = "****"
        elseif pval <= 0.001
            tag = "***"
        elseif pval <= 0.01
            tag = "**"
        elseif pval <= 0.05
            tag = "*"
        else
            tag = "ns"
        end
        
        print(pval)
   
        
        annotate(string("P = ",round(pval, sigdigits = 2
                    )), xy=(mean([0, xs[i]]), ys[i]+maximum(pd[string(col)])*0.05),
        horizontalalignment="center",
        verticalalignment="center")
        
        annotate(100 - round(mean(ko)*100, digits = 2), xy=(xs[i], 10),
        horizontalalignment="center",
        verticalalignment="center")
    
   

        pretty_axes()
        
    end
    
    ylabel("% of Rad21 deleted", fontsize = 14)
    
    plt.tight_layout()

xticks([0, 1, 2, 3], ["-4OHT \n 3 days", "+4OHT \n 1 day", "+4OHT \n 2 days", "+4OHT \n 3 days"], fontsize = 11)
    yticks(fontsize = 11)

xlabel("Rad21 lox/lox \n macrophages", fontsize = 14)
    pretty_axes2()

    
savefig("genetic_Rad21_left.png")
savefig("genetic_Rad21_left.svg")
end