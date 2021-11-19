

function make_summary_bursts_lps(GENE;limit = 1, rep = 1)
    set = sort!(get_genedata2(GENE), :Timepoint)
    
    
    
    reps = split_by(set, :Rep)
    
    newdf =  vcat([calculate_bf(reps[ii], :Sample; limit = limit, rep = ii) for ii in keys(reps)]...)
    
    newdf[!,:Gene] =[GENE for ii in 1:nrow(newdf)]
    
    newdf[!,:Genotype] = [split(ii, "_")[1] for ii in newdf[!,:Sample]]
    
    newdf[!,:Timepoint] = [split(ii, "_")[2] for ii in newdf[!,:Sample]]
    
    newdf = sort(newdf, :Rep, rev = true)
     newdf = sort(newdf, :Genotype, rev = true)
     newdf = sort(newdf, :Timepoint, rev = true)

end

function get_genedata2(genesymbol)
    genedata = CSV.read(ENV["Code"]* "/../Code_Paper/CompleteSets/CompleteSets/"*genesymbol*".csv", DataFrames.DataFrame)
       try
        genedata[!,:TSS1_r2] = [if ii == "NA" 0 else parse(Float64, ii) end for ii in genedata[!,:TSS1_r2]]
    catch
    end
    
       try
        genedata[!,:TSS2_r2] = [if ii == "NA" 0 else parse(Float64, ii) end for ii in genedata[!,:TSS2_r2]]
    catch
    end
    genedata[!,:Sample] = genedata[!,:Genotype] .* "_" .* string.(genedata[!,:Timepoint])
    genedata
end



function calculate_bf(df, col; limit = 1, rep = 1)
   sams_df = split_by(df, col)
    new_df = DataFrames.DataFrame(
        Sample = collect(keys(sams_df)),
        N_Cells = [nrow(sams_df[ii]) for ii in keys(sams_df)],
        N_TSS =  [nrow(tss_data(sams_df[ii], limit = limit)) for ii in keys(sams_df)],
        Mean_TSS =  [mean(tss_data(sams_df[ii], limit = limit)[!,:TSS1_r2]) for ii in keys(sams_df)],
        Median_TSS =  [try median(tss_data(sams_df[ii], limit = limit)[!,:TSS1_r2]) catch; NaN end for ii in keys(sams_df)]
    )
    
    new_df[!,:BF] = new_df[!,:N_TSS]./new_df[!,:N_Cells]./2
    new_df[!,:Rep] = [rep for ii in 1:nrow(new_df)]
    new_df
    
end
                                                          
limit = 1
a = make_summary_bursts_lps("Peli1", limit = limit)
a[!,:Timepoint] = [if ii .== "135" "120" else ii end for ii in a[!,:Timepoint]]
a[!,:Sample] = a[!,:Genotype] .* "_" .* a[!,:Timepoint]

all_genes = vcat(
    make_summary_bursts_lps("Ifnb1", limit = limit),
    make_summary_bursts_lps("Il12b", limit = limit),
    make_summary_bursts_lps("Ifit1", limit = limit),
    make_summary_bursts_lps("Cxcl10", limit = limit),
    make_summary_bursts_lps("Egr2", limit = limit),
    make_summary_bursts_lps("Prdm1", limit = 0),
    make_summary_bursts_lps("Sertad2", limit = limit),
   a,
    make_summary_bursts_lps("Fh1", limit = limit),
    make_summary_bursts_lps("Hprt", limit = limit),
)


CSV.write("all_genes.csv", all_genes)

