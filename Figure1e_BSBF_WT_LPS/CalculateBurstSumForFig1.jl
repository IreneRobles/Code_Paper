function calculate_bf(t, tname; limit = 0)
    new_df = DataFrames.DataFrame()
    
    samples = unique(t[!,:Sample])
    new_df[!,:Sample] = samples
    new_df[!,:Genotype] = [split(ii, "_")[1] for ii in samples]
    new_df[!,:Timepoint] = [split(ii, "_")[2] for ii in samples]
    
    n_cells = []
    cellsactive1 = []
    cellsactive2 = []
    n_tss =[]
    bs_mean =[]
    bs_median =[]
     bs_std =[]
    
    for ii in samples
        push!(n_cells, sum(t[!,:Sample] .== ii))
        sp_sam = t[t[!,:Sample] .== ii, :]
        locus1_bool = sp_sam[!,Symbol(string("TSS1_r2"))] .> limit
        locus1 = sp_sam[locus1_bool,Symbol(string("TSS1_r2"))]
        locus2_bool = sp_sam[!,Symbol(string("TSS2_r2"))] .> limit
        locus2 = sp_sam[locus2_bool,Symbol(string("TSS2_r2"))]
        push!(cellsactive1, Statistics.mean(locus1_bool.|locus2_bool))
        push!(cellsactive2, Statistics.mean(locus1_bool.*locus2_bool))
        push!(n_tss, length(locus1)+length(locus2))
        sizes = append!(locus1, locus2)
        push!(bs_mean, Statistics.mean(sizes))
        mediana = if isempty(sizes) 0 else median(sizes) end
        push!(bs_median, mediana)
        push!(bs_std, Statistics.std(sizes))
        
    end
    
    new_df[!,:N_cells] = n_cells
    new_df[!,Symbol(string("CellsActive1_", tname))] = cellsactive1
    new_df[!,Symbol(string("CellsActive2_", tname))] = cellsactive2
    new_df[!,Symbol(string(tname, "_N"))] = n_tss
    
    new_df[!,Symbol(string("BF_", tname))] = n_tss./2n_cells
    new_df[!,Symbol(string("BS_mean_", tname))] = bs_mean
    new_df[!,Symbol(string("BS_median_", tname))] = bs_median
    new_df[!,Symbol(string("BS_std_", tname))] = bs_std
    
    new_df
    
end
      

BF = Dict()
bf = calculate_bf_by_rep(get_completeset("Peli1"), "TSS", "Peli1")
# Only consider 1 timepoint
bf[!,:Timepoint] = [if ii .== "135" "120" else ii end for ii in bf[!,:Timepoint]]
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "120"
bool = bool1.| bool2
BF["Peli1"] = bf[bool, :]

bf = calculate_bf_by_rep(get_completeset("Il12b"), "TSS", "Il12b")
# Only consider 1 timepoint
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "90"
bool = bool1.| bool2
BF["Il12b"] = bf[bool, :]


bf = calculate_bf_by_rep(get_completeset("Ifnb1"), "TSS", "Ifnb1")
# Only consider 1 timepoint
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "90"
bool = bool1.| bool2
BF["Ifnb1"] = bf[bool, :]


bf = calculate_bf_by_rep(get_completeset("Ifit1"), "TSS", "Ifit1")
# Only consider 1 timepoint
BF["Ifit1"] = bf

bf = calculate_bf_by_rep(get_completeset("Cxcl10"), "TSS", "Cxcl10")
# Only consider 1 timepoint
BF["Cxcl10"] = bf
            
bf = calculate_bf_by_rep(get_completeset("Sertad2"), "TSS", "Sertad2")
# Only consider 1 timepoint
BF["Sertad2"] = bf

bf = calculate_bf_by_rep(get_completeset("Prdm1"), "TSS", "Prdm1")
# Only consider 1 timepoint
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "60"
bool = bool1.| bool2
BF["Prdm1"] = bf[bool, :]

bf = calculate_bf_by_rep(get_completeset("Egr2"), "TSS", "Egr2")
# Only consider 1 timepoint
BF["Egr2"] = bf


bf = calculate_bf_by_rep(get_completeset("Fh1"), "TSS", "Fh1")
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "8"
bool = bool1.| bool2
# Only consider 1 timepoint
BF["Fh1"] = bf[bool, :]


bf = calculate_bf_by_rep(get_completeset("Hprt"), "TSS", "Hprt")
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "8"
bool = bool1.| bool2
# Only consider 1 timepoint
BF["Hprt"] = bf[bool, :]


BFs = join_in_all_common_columns([BF[key] for key in keys(BF)]...)

CSV.write("BurstFrequencies.csv", BFs)
                                                            
  