ENV["Code"] = "../../Code"
for folder in readdir(ENV["Code"]); push!(LOAD_PATH, normpath(ENV["Code"], folder));end
include(ENV["Code"]*"/../Code_Paper/Code/meanmRNAcounts_BSBF.jl")
using Random, Statistics
using RCall


function bootstrap_mean_std(genedata, gene; subsampling = true, cellssubsampled = 50, times = 1000, limit = 2, limit_exon = 5)
    col = :Sample_Rep
    if !in( "N_exon",names(genedata))
        genedata[!,"N_exon"] = genedata[!,:TSS1_r2].+ genedata[!,:TSS2_r2] .+ genedata[!,:N_thres_Total]
    end
    new_df = DataFrames.DataFrame()
    genedata[!,:Sample] = genedata[!,:Genotype] .* "_" .* string.(genedata[!,:Timepoint])
    genedata[!,:Sample_Rep] = genedata[!,:Sample] .* "_" .* string.(genedata[!,:Rep])
    samples = unique(genedata[!,col])
    new_df[!,:Sample_Rep] = samples
    new_df[!,:Genotype] = [split(ii, "_")[1] for ii in samples]
    new_df[!,:Timepoint] = [split(ii, "_")[2] for ii in samples]
    new_df[!,:Rep] = [split(ii, "_")[3] for ii in samples]
    new_df[!,:Sample] = new_df[!,:Genotype] .* "_" .* new_df[!,:Timepoint]
    
    spt = split_by(genedata, col)
    means = []
    stds = []
    means_ = []
    stds_ = []
    bs = []
    bf = []
    ns = []
    nbursts = []
    exon_fraction_cells = []
    exon_pos_cells = []
    tss_fraction_cells = []
    tss_pos_cells = []
    
    
    for ii in samples
        distribution = spt[ii][!,:N_exon]
        n = length(distribution)
        push!(ns, n)
        subdist = []
        if subsampling
            subdist  = [distribution[shuffle(1:n)[1:cellssubsampled]] for jj in 1:times]
        else
            subdist = distribution
        end
        push!(means, Statistics.mean([Statistics.mean(jj) for jj in subdist]))
        push!(stds, Statistics.mean([Statistics.std(jj) for jj in subdist]))
        push!(means_, Statistics.mean(distribution))
        push!(exon_fraction_cells, Statistics.mean(distribution .>= limit_exon))
        push!(exon_pos_cells, Statistics.sum(distribution .>= limit_exon))
        push!(stds_, Statistics.std(distribution))
        distribution = vcat(spt[ii][!,:TSS1_r2], spt[ii][!,:TSS2_r2])
        push!(bs, Statistics.mean(distribution[distribution.>limit]))
        push!(bf, Statistics.mean(distribution.>limit))
        push!(nbursts, Statistics.sum(distribution.>limit))
        bool1 = spt[ii][!,:TSS1_r2] .> limit
        bool2 = spt[ii][!,:TSS2_r2] .> limit
        bool = bool1 .+ bool2
        push!(tss_pos_cells, Statistics.sum(bool.>0 ))
        push!(tss_fraction_cells, Statistics.mean(bool.>0 ))

    end
    new_df[!,:n_cells] = ns
    new_df[!,:nburst] = nbursts
    new_df[!,:Gene] = [gene for ii in means_]
    new_df[!,:mean] = means_
    new_df[!,:bootstrap_mean] = means
    new_df[!,:std] = stds_
    new_df[!,:bootstrap_std] = stds
    new_df[!,:fractionTSS] = tss_fraction_cells
    new_df[!,:fractionExon] = exon_fraction_cells
    new_df[!,:posTSS] = tss_pos_cells
    new_df[!,:posExon] = exon_pos_cells
    new_df[!,"TSS_BS"] =  bs
    new_df[!,"momment_BS"] =  (new_df[!,"bootstrap_std"].^2)./ new_df[!,"bootstrap_mean"]
    new_df[!,"TSS_BF"] =  bf
    new_df[!,"momment_BF"] = new_df[!,"bootstrap_mean"] ./ ( new_df[!,"bootstrap_std"].-1)
    
    new_df
    
end