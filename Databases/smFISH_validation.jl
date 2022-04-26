module smFISH_validation

folder = "../Databases"

using CSV
using DataFrames

function smFISHdata(gene)
    GENE = uppercase(gene)
    if in(GENE, ["ACTB", "IL12B", "STAT1"])
        return CSV.read(smFISH_validation.folder*"/smFISH_validation/actb_stat1_il12b.csv", DataFrames.DataFrame)
    elseif in(GENE, ["ATF3", "TNF", "CD40"])
        return CSV.read(smFISH_validation.folder*"/smFISH_validation/cd40_atf3_tnf.csv", DataFrames.DataFrame)
    end
end

function fano(vector)
    std(vector)^2/mean(vector)
end

function variation(vector)
     std(vector)/mean(vector)
end

function cv2(vector)
    std(vector)^2/mean(vector)^2
end

function apply_func_per_sample_replicate(tb, funct, functname)
    reps = unique(tb[!,:Rep])
    genes = names(tb)[6:8]
    rep_sam =  Array{String,1}()
    gene_col = Array{String,1}()
    datapoint_col = Array{Float64,1}()
    
    for rep in reps
        tb_rep = tb[tb[!,:Rep].== rep, :]
        samples = unique(tb_rep[!,:Sample])
        for sam in samples
            tb_sam = tb_rep[tb_rep[!,:Sample].== sam, :]
            
            for gene in genes
                push!(rep_sam, string(rep, "__", sam))
                push!(gene_col, string(gene))
                push!(datapoint_col, funct(tb_sam[:,gene]))
                
            end
        end
    end
    df = DataFrame()
    df[!,:Rep_Sample] = rep_sam
    df[!,:Rep] = [split(ii, "__")[1] for ii in df[!,:Rep_Sample]] 
    df[!,:Sample] = [split(ii, "__")[2] for ii in df[!,:Rep_Sample]] 
    df[!,:Genotype] = [split(ii, "_")[3] for ii in df[!,:Rep_Sample]] 
    df[!,:Timepoint] = [split(ii, "_")[4] for ii in df[!,:Rep_Sample]] 
    df[!,:Gene] = gene_col
    df[!,Symbol(functname)] = datapoint_col
    
    df
end

    


end