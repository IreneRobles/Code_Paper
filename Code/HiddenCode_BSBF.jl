using Seaborn, NoLongerProblems_Pandas, Pandas, NoLongerProblems, PrettyPlotting, DataFrames, PyPlot
using NoLongerProblems, Statistics
using CSV, DataFrames, MultipleTesting
using FQfiles, RCall, PyCall
    

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
    
    try
        new_df[!,:Timepoint] = map(x-> parse(Int, split(x, "_")[2]), new_df[!,:Sample])
    catch
        new_df[!,:Timepoint] = map(x-> split(x, "_")[2], new_df[!,:Sample])
    end

    
 function calculate_burst_size(sample; func = mean)
        d = completecases!(df[df[!,:Sample].==sample, :])
        try
        burst_data = collect(vcat(d[!,Symbol(string("TSS1_", int))], d[!,Symbol(string("TSS2_", int))]));
            catch; 
            burst_data = collect(vcat(d[!,Symbol(string("locus1_", int))], d[!,Symbol(string("locus2_", int))]));
        end
            
         bool = burst_data .> limit
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
        bool = burst_data .> limit
        burst_data = length(burst_data[bool])/n
    end
    
    function fraction_expressing(sample, limit)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)
        
        o = d[!,:TSS1_r2].+ d[!,:TSS2_r2] .> limit
        
        p = d[!,:N_total] .> limit
     
        burst_data = p .+  o;
  
        bool = burst_data .> 0
        burst_data = length(burst_data[bool])/n
    end
    
        
    function fraction_bursting(sample, limit)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)
        
        o = d[!,:TSS1_r2].+ d[!,:TSS2_r2] 
        
     
        burst_data =  o;
  
        bool = burst_data .> limit
        burst_data = length(burst_data[bool])/n
    end
    
        function fraction_expressing_n(sample, limit)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)
        
        o = d[!,:TSS1_r2].+ d[!,:TSS1_r2] .> limit
        
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
  
        bool = burst_data .> limit
        burst_data = length(burst_data[bool])
    end
    
        
    function fraction_bursting_n(sample, limit)
        d = completecases!(df[df[!,:Sample].==sample, :])
        n = nrow(d)
        
        o = d[!,:TSS1_r2].+ d[!,:TSS1_r2] .> limit
        
     
        burst_data =  o;
  
        bool = burst_data .> limit
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
        bool = burst_data .> limit
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


py"""
import matplotlib as mpl

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
mpl.rcParams.update(new_rc_params)
import matplotlib.pyplot as plt
"""
function get_sample(df, sample)
    f(x) = [i == sample for i in df[:Sample]]
    df[f(df), :]
end

function positive_cells_(df, gen, timepoint, thing; return_ = "trues", limit = 1)
    data_to_plot = get_sample(df, string(gen, "_", timepoint))
    
     data_to_plot = data_to_plot[thing]
    
    total = length(data_to_plot)
    
     bool = [i .>= limit for i in data_to_plot]
    
    trues = sum(bool)
    f_trues = trues / total
    
    if return_ == "trues"
        return f_trues
    else
        return 1-f_trues
    end
    

end

function CV(df, gen, timepoint, thing; return_ = "trues", limit = 1)
    data_to_plot = get_sample(df, string(gen, "_", timepoint))
    
     data_to_plot = data_to_plot[thing]
    
    total = length(data_to_plot)
    
     bool = std(data_to_plot)/mean(data_to_plot)

end

function Mean(df, gen, timepoint, thing; return_ = "trues", limit = 1)
    data_to_plot = get_sample(df, string(gen, "_", timepoint))
    
     data_to_plot = data_to_plot[thing]
    
    total = length(data_to_plot)
    
     bool = mean(data_to_plot)

end

function Median(df, gen, timepoint, thing; return_ = "trues", limit = 1)
    data_to_plot = get_sample(df, string(gen, "_", timepoint))
    
     data_to_plot = data_to_plot[thing]
    
    total = length(data_to_plot)
    
     bool = median(data_to_plot)

end


function plot_fraction_of_positive_cells(df, key; timepoints = [0, 60, 90, 120], color = ["blue", "orange"], rep = 1, limit = 2)

    wt = [positive_cells_(df, "WT", i, key, limit = limit) for i in timepoints]
    ko = [positive_cells_(df, "Rad21KO", i, key, limit = limit) for i in timepoints]
    PyPlot.scatter(timepoints, wt, c = color[1])
    PyPlot.scatter(timepoints, ko, c = color[2])
    
    PyPlot.plot(timepoints, wt, label = string("WT_", rep), c = color[1])
    PyPlot.plot(timepoints, ko, label = string("Rad12KO_", rep), c = color[2])

    pretty_axes()
end

function fraction_of_positive_cells(df::DataFrames.DataFrame, key; timepoints = [0, 60, 90, 120], color = ["blue", "orange"], rep = 1, limit = 2)
    new_df = DataFrames.DataFrame()
    new_df[!,:Sample] = unique(df[:Sample])
    new_df[!,:Genotype]  = [split(i, "_")[1] for i in new_df[:Sample]]
    new_df[!,:Timepoint]  = [parse(Int, split(i, "_")[2]) for i in new_df[:Sample]]
    new_df[!,:F_Pos]  = [positive_cells_(df, new_df[i, :Genotype], new_df[i, :Timepoint], key, limit = limit) for i in 1:nrow(new_df)]
    new_df[!,:CV]  = [CV(df, new_df[i, :Genotype], new_df[i, :Timepoint], key, limit = limit) for i in 1:nrow(new_df)]
    new_df[!,:Mean]  = [Mean(df, new_df[i, :Genotype], new_df[i, :Timepoint], key, limit = limit) for i in 1:nrow(new_df)]
    new_df[!,:Median]  = [Median(df, new_df[i, :Genotype], new_df[i, :Timepoint], key, limit = limit) for i in 1:nrow(new_df)]

    new_df[!,:Replicate] = [rep for i in  1:nrow(new_df)]
    
    new_df
   
end

function fraction_of_positive_cells(dfs::Array, key; timepoints = [0, 60, 90, 120], color = ["blue", "orange"], limit = 2)
    i = 0
    ds = []
    for df in dfs
        i+=1
        d = fraction_of_positive_cells(df, key; timepoints = timepoints, rep = i, limit = limit)
        push!(ds, d)
    end
    
    ds = join_in_all_common_columns(ds...)
    
    pd = Pandas.DataFrame(ds)
    
    Seaborn.boxplot(data = pd, x = "Timepoint", y = "F_Pos", hue = "Genotype")
    ylabel("Fraction of cells with \n more than $limit counts", fontsize = 20)

end






function do_mantelhaen(df, comp1, comp2)
    s1 = df[df[!,:Sample] .== comp1, :]
    s2 = df[df[!,:Sample] .== comp2, :]
    tb_test = join_in_all_common_columns(s1, s2)
    tb1 = tb_test; tb1[!,:Burst] = ["Yes" for ii in 1:nrow(tb_test)]
    tb1[!,:Count] = tb1[!,:nburst]
    tb1 = tb1[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb2 = tb_test; tb2[!,:Burst] = ["No" for ii in 1:nrow(tb_test)]
    tb2[!,:Count] = tb2[!,:n_cells] * 2 - tb2[!,:nburst]
    tb2 = tb2[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb = join_in_all_common_columns(tb1, tb2)

    
test = R"""
library("psych")
library("vcd")
library("DescTools")
library("rcompanion")
library("dplyr")
options(warn=-1)

tb = $tb
tb$Count = as.integer(tb$Count)
    
Data <- mutate(tb,
           Sample = factor(Sample, levels=unique(Sample)),
           Burst = factor(Burst, levels=unique(Burst)),
           Rep = factor(Rep, levels=unique(Rep))
           )

# Last variable is the strata (the variable that is not check for assotiation)
Data.xtabs <- xtabs(Count ~ Burst + Sample + Rep, 
                       data=Data)

ftable(Data.xtabs)  

mantelhaen.test(Data.xtabs)
"""

return test[3][1],test[5][1]
end




