ENV["Code"] = "../../Code"
for folder in readdir(ENV["Code"]); push!(LOAD_PATH, normpath(ENV["Code"], folder));end

using CSV, DataFrames
using NoLongerProblems_FileHandling, NoLongerProblems, NoLongerProblems_Pandas
using MultipleTesting, HypothesisTests, Statistics
using Seaborn, PyPlot, PyCall, RCall
import Pandas
using Random, PrettyPlotting


function tss_data(tss; limit = 1)
    tss1 = tss[tss[!,"TSS1_r2"].>limit, :]
    tss2 = tss[tss[!,"TSS2_r2"].>limit, :]
    tss2[!,"TSS1_r2"] = tss2[!,"TSS2_r2"]
    tsss = vcat(tss1, tss2)
end

function get_completeset(genename)
    tb = CSV.read(string(ENV["Code"], "/../Code_Paper/CompleteSets/CompleteSets/",genename,".csv"),DataFrames.DataFrame)
    tb[!,:Gene]= [genename for ii in 1:nrow(tb)]
    tb = sort!(tb, :Genotype, rev = true)
    tb = sort!(tb, :Timepoint, rev = false)

    return tb
end


function YiFangDeseq(name)
    yifandeseqfolder = ENV["Code"]*"/Databases/Cuartero2018/DE_result_norm2Rep_n_Spikes/"
    files = readdir(yifandeseqfolder)
     bool = [occursin(name, ii) for ii in files]
    file = files[bool]
     if length(file) == 1
        file = yifandeseqfolder*file[1]
        d = dropmissing(CSV.read(file, DataFrames.DataFrame, missingstring = "NA"))
        return rename!(d, :ensembl_id => :EnsemblID, :mgi_symbol => :GeneSymbol)
    else
        return files
    end
end

function CSV.read(file, typeData;kwards...)
    typeData(CSV.File(file; kwards...))
end

function GroseDeseq(name)
    folder = "../../MiniPrograms_Paper1/GROseq_DEseq/"

    files = readdir(folder)
     bool = [occursin(name, ii) for ii in files]
    file = files[bool]
     if length(file) == 1
        file = folder*file[1]
         d = innerjoin(dropmissing(CSV.read(file, DataFrames.DataFrame, missingstring = "NA")), mm9.EnsemblIDtoGeneSymbol_table(), on = :EnsemblID)
        return  d
    else
        return ""
    end
end

function dysplay_correlation(tb, x, y; style = "RNAseqcomp", hue = "GeneSymbol", palette = nothing, xy = (0, 4))
    bool1 = .! isnan.(tb[!,x]); bool2 = .! isnan.(tb[!,y])
    nonantb = tb[bool1.*bool2, :]
    sort!(nonantb, "smFISHcomp")
    x1 = [ii for ii in nonantb[!,x]]
    y1 = [ii for ii in nonantb[!,y]]
    
    pdt = Pandas.DataFrame(sort(nonantb, :Genotype, rev = true))

    
    py"""
    import seaborn as sns
    sns.scatterplot(data = $pdt, x= $x, y = $y, hue = $hue, style = $style, palette = $palette, s = 50, linewidth = 0)
    """
    
    pretty_axes2()
    legend_out_of_plot()
    test = R"cor.test($x1,$y1)"
    pval = round(test[3][1], sigdigits = 2)
    corr = round(test[4][1], sigdigits = 2)
    annotate("""
        r = $corr
        pval = $pval
        """, xy = xy, va = "center"
    )
    
    
end

function calculate_bf(t, tname; limit = 2)
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
    exon_mean = []
    
    for ii in samples
        push!(n_cells, sum(t[!,:Sample] .== ii))
        sp_sam = t[t[!,:Sample] .== ii, :]
        try 
            push!(exon_mean, Statistics.mean(sp_sam[!,"N_exon"]))
        catch
            push!(exon_mean, Statistics.mean(sp_sam[!,"N_thres_Total"]))
        end
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
    new_df[!,Symbol(string("log2BF_", tname))] = log2.(n_tss./2n_cells)
    new_df[!,Symbol(string("BS_mean_", tname))] = bs_mean
    new_df[!,Symbol(string("BS_median_", tname))] = bs_median
    new_df[!,Symbol(string("BS_std_", tname))] = bs_std
    new_df[!,Symbol(string("MeanCounts", tname))] = exon_mean
                    
    
    new_df
    
end
function calculate_bf_by_rep(df, tname, genename; limit = 1)
    
    df[!,:Sample] = [string(df[ii, :Genotype], "_",df[ii, :Timepoint]) for ii in 1:nrow(df)]
    subdfs = NoLongerProblems.split_by(df, :Rep)
    reps = collect(keys(subdfs))
    
    subdfs = [calculate_bf(subdfs[ii], tname; limit = limit) for ii in reps]
    
    
    for ii in 1:length(reps)
          subdfs[ii][!,:Rep]= [reps[ii] for a in 1:nrow(subdfs[ii])]
        subdfs[ii][!,:Gene]= [genename for a in 1:nrow(subdfs[ii])]
    end
    
    tb = join_in_all_common_columns(subdfs)
                    
    CSV.write(ENV["Code"]*"/../Code_Paper/CompleteSets/"*"BurstSummary/"*genename*".csv", tb)
                    
    return tb
end

function calculate_bf_by_rep(df, tname, genename; limit = 2)
    
    df[!,:Sample] = [string(df[ii, :Genotype], "_",df[ii, :Timepoint]) for ii in 1:nrow(df)]
    subdfs = NoLongerProblems.split_by(df, :Rep)
    reps = collect(keys(subdfs))
    
    subdfs = [calculate_bf(subdfs[ii], tname; limit = limit) for ii in reps]
    
    
    for ii in 1:length(reps)
          subdfs[ii][!,:Rep]= [reps[ii] for a in 1:nrow(subdfs[ii])]
        subdfs[ii][!,:Gene]= [genename for a in 1:nrow(subdfs[ii])]
    end
    
    join_in_all_common_columns(subdfs)
end

function plot_bf(gene; scale = "min")

pd = Pandas.DataFrame(BF[gene])
Seaborn.boxplot(data = pd, y = "BF_TSS", hue = "Genotype", x = "Timepoint")
Seaborn.stripplot(data = pd, y = "BF_TSS", hue = "Genotype", x = "Timepoint", dodge = true)
pretty_axes2()
xlabel("Time after LPS ("*scale*")")
ylabel("Burst Freq.")
    title(gene)
    legend_removal()
    
end

function plot_bs(gene, limit, time1, time2; scale = "min", maxy = 50)
df = tss_data(get_completeset(gene); limit = limit)
df[!,:Timepoint] = [if ii .== "135" "120" else ii end for ii in df[!,:Timepoint]]
bool1 = df[!,:Timepoint] .== time1 
bool2 = df[!,:Timepoint] .== time2
bool = bool1 .| bool2
pd = Pandas.DataFrame(df[bool, :])
Seaborn.boxplot(data = pd, y = "TSS1_r2", hue = "Genotype", x = "Timepoint", showfliers = false)
Seaborn.stripplot(data = pd, y = "TSS1_r2", hue = "Genotype", x = "Timepoint", dodge = true, zorder = 0)
pretty_axes2()
xlabel("Time after LPS ("*scale*")")
ylabel("Burst Size")
    title(gene)
    legend_removal()
    ylim(0, maxy)
    
end


