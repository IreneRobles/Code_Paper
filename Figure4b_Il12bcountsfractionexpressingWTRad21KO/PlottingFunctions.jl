function Il12b_maturefig_expressing()

    pd = Pandas.DataFrame(reps_expressingexon)
     Seaborn.boxplot(data = pd, x= "Timepoint", y = "N_exon", palette = ["darkgray", "red"], hue = "Genotype", showfliers = false, )
    #Seaborn.stripplot(data = pd, x= "Timepoint", y = "N_exon", hue = "Rep",palette = "Greys", size = 1.5,jitter = 0.45, zorder = 0)
    pretty_axes2()
    ylim(0, 150)

    h = 35
    u = 7

    xy = [-0.25, 0.25]
    plt.plot(xy, [h, h], c = "black", lw = 0.75)
    plt.annotate(string("", round(t[t[!,"Column1"].== "WT_0-Rad21KO_0", "p adj"][1], sigdigits = 2)), xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    h = 70
    xy = [0.75, 1.25]
    plt.plot(xy, [h, h], c = "black", lw = 0.75)
    plt.annotate(string("", round(t[t[!,"Column1"].== "WT_60-Rad21KO_60", "p adj"][1], sigdigits = 2)), xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    h = 123
    xy = [1.75, 2.25]
    plt.plot(xy, [h, h], c = "black", lw = 0.75)
    plt.annotate(string("", round(t[t[!,"Column1"].== "WT_90-Rad21KO_90", "p adj"][1], sigdigits = 2)), xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    h = 130
    xy = [2.75, 3.25]
    plt.plot(xy, [h, h], c = "black", lw = 0.75)
    plt.annotate(string("", round(t[t[!,"Column1"].== "WT_120-Rad21KO_120", "p adj"][1], sigdigits = 2)), xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    
    
    ylabel("Il12b  mature mRNA counts \n in cells with more than $lim_exp")
    xlabel("Time after LPS (min)")
    
    squareplot()
    line075black()
    
    
end

function locusdata2(tab...; genotype = :Genotype, timepoint = true, 
        tss1 = :TSS1_r2, 
        tss2 = :TSS2_r2
    )
    new_df = DataFrames.DataFrame()
    for ii in 1:length(tab)
           if timepoint == false
                return tab[!,ii][!,:Timepoint]  = [0 for ii in 1:nrow(tab[!,ii])]
            end
         tab[ii][!,:Genotype]  = tab[ii][!,genotype]
        
        tss1s =  tab[ii][!, [:Cell, :Sample, :Genotype, :Timepoint, :Image, tss1]]; 
        rename!(tss1s, tss1 => :BurstSize)
        tss2s =  tab[ii][!, [:Cell, :Sample, :Genotype, :Timepoint, :Image, tss2]]; 
        rename!(tss2s, tss2 => :BurstSize)
        
        new_tb = join_in_all_common_columns(tss1s, tss2s)
        new_tb[!,:Rep] = [ii for a in 1:nrow(new_tb)]
        if ii == 1
            new_df = new_tb
        else
            new_df = join_in_all_common_columns(new_df, new_tb)
        end
        
    end
    return new_df
end

completecases! = dropmissing!

function mean_burst_size_and_burst_fraction(df; 
        int = :integrated_intensity_ratio_With_radious_1_around_brightest_pixel_, 
        rep = 1,
        limit = 2
    )
    
    #df[:Total_int] =  df[Symbol(string("TSSsum_", int))] + df[:N_intron_Nuc]
    
    new_df = DataFrames.DataFrame()
    new_df[!,:Sample] = unique(df[!,:Sample])
    new_df[!,:Genotype] = map(x-> split(x, "_")[1], new_df[!,:Sample])
      new_df[!,:Timepoint] = map(x-> parse(Int, split(x, "_")[2]), new_df[!,:Sample])
    
 function calculate_burst_size(sample; func = mean)
        d = completecases!(df[df[!,:Sample].==sample, :])
        try
        burst_data = collect(vcat(d[!,Symbol(string("TSS1_", int))], d[!,Symbol(string("TSS2_", int))]));
            catch; 
            burst_data = collect(vcat(d[!,Symbol(string("locus1_", int))], d[!,Symbol(string("locus2_", int))]));
        end
            
         bool = burst_data .!= 0.0
        burst_data =  burst_data[bool]
        try func(burst_data) catch; 0 end
    end

     function calculate_burst_fraction(sample)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)*2
       try
        burst_data = collect(vcat(d[!,Symbol(string("TSS1_", int))], d[!,Symbol(string("TSS2_", int))]));
            catch;
            burst_data = collect(vcat(d[!,Symbol(string("locus1_", int))], d[!,Symbol(string("locus2_", int))]));
        end
        bool = burst_data .!= 0.0
        burst_data = length(burst_data[bool])/n
    end
    
    function fraction_expressing(sample, limit)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)
        
        o = d[!,:TSS1_r2].+ d[!,:TSS1_r2] .> 0
        
        p = d[!,:N_total] .> limit
     
        burst_data = p .+  o;
  
        bool = burst_data .> 0
        burst_data = length(burst_data[bool])/n
    end
    
        
    function fraction_bursting(sample, limit)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)
        
        o = d[!,:TSS1_r2].+ d[!,:TSS1_r2] .> 0
        
     
        burst_data =  o;
  
        bool = burst_data .> 0
        burst_data = length(burst_data[bool])/n
    end
    
        function fraction_expressing_n(sample, limit)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)
        
        o = d[!,:TSS1_r2].+ d[!,:TSS1_r2] .> 0
        
        p = d[!,:N_intron] .> limit
     
        burst_data = p .+  o;
  
        bool = burst_data .> 0
        burst_data = length(burst_data[bool])
    end
    
    function fraction_expressing_exon(sample, limit)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)
                
        p = d[!,:N_exon] .> limit
     
        burst_data = p 
  
        bool = burst_data .> 0
        burst_data = length(burst_data[bool])
    end
    
        
    function fraction_bursting_n(sample, limit)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)
        
        o = d[!,:TSS1_r2].+ d[!,:TSS1_r2] .> 0
        
     
        burst_data =  o;
  
        bool = burst_data .> 0
        burst_data = length(burst_data[bool])
    end

    function calculate_noburst(sample)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)*2
        try
        burst_data = collect(vcat(d[!,Symbol(string("TSS1_", int))], d[!,Symbol(string("TSS2_", int))]));
            catch 
            burst_data = collect(vcat(d[!,Symbol(string("locus1_", int))], d[!,Symbol(string("locus2_", int))]));
        end
        bool = burst_data .!= 0.0
        burst_data = n - length(burst_data[bool])
    end
    
     function f(sample, col, func)
        d = completecases!(df[df[!,:Sample].==sample, :])
        return func(d[!,col])
    end
    
     function countmore(sample, col, limit)
        d = completecases!(df[df[!,:Sample].==sample, :])
        return count( x -> x> limit , d[!,col])/nrow(d)
    end
    
    
     new_df[!,:BurstSize_avg] = map(x-> calculate_burst_size(x, func = mean), new_df[!,:Sample])
    new_df[!,:nburst] = map(x-> calculate_burst_size(x, func = length), new_df[!,:Sample])
     new_df[!,:nburst_no] = map(x-> calculate_noburst(x), new_df[!,:Sample])
    new_df[!,:BurstSize_median] = map(x-> calculate_burst_size(x, func = median), new_df[!,:Sample])
    new_df[!,:BurstFraction] = map(x-> calculate_burst_fraction(x), new_df[!,:Sample])
    
    new_df[!,:N_Expressing] = map(x-> fraction_expressing_n(x, limit), new_df[!,:Sample])
    new_df[!,:N_Expressing_exon] = map(x-> fraction_expressing_exon(x, limit), new_df[!,:Sample])

    new_df[!,:N_Bursting] = map(x-> fraction_bursting_n(x, limit), new_df[!,:Sample])
   
    
    
    new_df[!,:Expressing] = map(x-> fraction_expressing(x, limit), new_df[!,:Sample])
    new_df[!,:Bursting] = map(x-> fraction_bursting(x, limit), new_df[!,:Sample])
    
    new_df[!,:n_cells] = (new_df[!,:nburst] .+ new_df[!,:nburst_no])/2
    
    new_df[!,:Rep] = [string(rep) for ii in 1:nrow(new_df)]
    
     
        
    new_df
    
end


function mean_burst_size_and_burst_fraction(df...; 
        limit = 0,
        int = :r2, 
    )
    
    join_in_all_common_columns([mean_burst_size_and_burst_fraction(df[ii], int = int, rep = ii, limit = limit) for ii in 1:length(df)])
    
end

function Il12b_Fexpcells_radwt()

    pd = Pandas.DataFrame(summary1)
     Seaborn.boxplot(data = pd, x= "Timepoint", y = "F_Expressing_exon", palette = ["darkgray", "red"], hue = "Genotype", showfliers = false, )
   # Seaborn.stripplot(data = pd, x= "Timepoint", y = "BurstFraction", hue = "Rep",palette = "Greys", size = 5,jitter = 0.1, zorder = 1)
    pretty_axes2()


    h = 0.1
    u = 0.03

    xy = [-0.25, 0.25]
    plt.plot(xy, [h, h], c = "black", lw = 0.75)
    plt.annotate(string("", round(ps[1], sigdigits = 2)), xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    h = 0.2
    xy = [0.75, 1.25]
    plt.plot(xy, [h, h], c = "black", lw = 0.75)
    plt.annotate(string("", round(ps[2], sigdigits = 2)), xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    h = 0.4
    xy = [1.75, 2.25]
    plt.plot(xy, [h, h], c = "black", lw = 0.75)
    plt.annotate(string("", round(ps[3], sigdigits = 2)), xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    h = 0.6
    xy = [2.75, 3.25]
    plt.plot(xy, [h, h], c = "black", lw = 0.75)
    plt.annotate(string("", round(ps[4], sigdigits = 2)), xy = [Statistics.mean(xy), h + u], ha = "center", va = "center")

    ylabel("Fraction of cells with more than \n $lim_exp Il12b mature mRNA counts ")


    xlabel("Time after LPS (min)")
    
        squareplot()
    line075black()
end

