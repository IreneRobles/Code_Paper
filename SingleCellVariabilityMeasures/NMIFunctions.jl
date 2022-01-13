
function calculate_entropy_per_cell(sceexp, assay)
    expmatrix = sceexp.assays[assay]
    ncells = size(expmatrix)[2]
    sceexp.colData[!, string(assay*"_Entropy")] = [get_entropy(expmatrix[:,ii]; base = 2, number_of_bins = 10) for ii in 1:ncells]
    sceexp
end

function calculate_mutualinformation(sceexp, assay)
    expmatrix = sceexp.assays[assay]
    ncells = size(expmatrix)[2]
    newassay = zeros(ncells, ncells)
    
    for ii in 1:ncells
        for jj in 1:ncells
            if ii != jj && newassay[ii, jj] == 0.0 
                mi  = get_mutual_information(expmatrix[:,ii], expmatrix[:,jj], base = 2, number_of_bins = 10)
                newassay[ii, jj] = mi
                newassay[jj, ii] = mi
            end
        end
    end
    sceexp.assays[string(assay, "_MI")] = newassay
    sceexp
end

function calculate_normalised_mutual_information(sceexp, assay)
    h = sceexp.colData[!,assay*"_Entropy"]
    mi = sceexp.assays[assay*"_MI"]
    expmatrix = sceexp.assays[assay]
    ncells = size(expmatrix)[2]
    newassay = zeros(ncells, ncells)
    for ii in 1:ncells
        for jj in (ii+1):ncells
            if ii != jj && newassay[ii, jj] == 0.0 
                nmi = mi[ii, jj] / sqrt(h[ii]*h[jj])
                newassay[ii, jj] = nmi
                newassay[jj, ii] = nmi
            end
        end
    end
    sceexp.assays[string(assay, "_NMI")] = newassay
    sceexp
end

function NMI_table(sceexp; assay = "CPM",times= 100, cellssubsampled = 100)
    samples = unique(sceexp.colData[!,:Sample])
    dfs = []
    for s in samples
        subsce = SingleCellExperiment.get_cells_with_this_characteristics([s], :Sample, sceexp)
        subsce = calculate_entropy_per_cell(subsce, assay)
        subsce = calculate_mutualinformation(subsce, assay)
        subsce = calculate_normalised_mutual_information(subsce, assay)
        # we subsampled 100 cells from each time-point (or mouse) 100 times and calculated the median NMI across each within-timepoint sampled pair.  
        nmi_median = zeros(times)
        n = size(subsce.assays[assay])[2]
        for aa in 1:times
            cols = shuffle(1:n)[1:cellssubsampled]
            nmis = []
            for ii in 1:cellssubsampled
                for jj in (ii+1):cellssubsampled
                        iiw = cols[ii]
                        jjw = cols[jj]
                        push!(nmis, subsce.assays[assay*"_NMI"][iiw,jjw])
                    end
                end
            nmi_median[aa] =  StatsBase.median(nmis)
        end
        new_df = DataFrames.DataFrame()
        new_df[!,:Sample] = [s for ii in 1:times]
        new_df[!,:NMI] = nmi_median
        push!(dfs, new_df)
    end
    
    return join_in_all_common_columns(dfs...)
    
end

function NMI_figure(table; u = 0.005, h = maximum(table[!,:NMI])+u, pvalue_ = :option2)
    sg = [split(ii, "_") for ii in table[!,:Sample]]
    table[!,:Genotype] = [ii[1] for ii in sg]
    table[!,:Timepoint] = [ii[2] for ii in sg]
    
    Seaborn.boxplot(data = Pandas.DataFrame(table), y = "NMI", x = "Timepoint", hue = "Genotype", palette = ["darkgray", "red"], showfliers = false)
    
    
    tims = unique(table[!,:Timepoint])
    # Do test
    for tim in 1:length(tims)
        plt.plot([tim-1.25, tim-0.75], [h,h])
        subtb = table[table[!,:Timepoint].==tims[tim], :]
        wt = subtb[subtb[!,:Genotype].=="WT", :NMI]
        rad = subtb[subtb[!,:Genotype].=="RAD21", :NMI]
        
        p = 1 
        
        if pvalue_ == :option1
            p = (sum(wt .< rad)+1)/((nrow(subtb)/2)+1)
        elseif pvalue_ == :option2
            p  = pvalue(MannWhitneyUTest(wt,rad))
        end
        println(p)
        x = (tim-1.25+tim-0.75)/2
        plt.annotate(round(p, sigdigits=2), xy = [x, h+u], va = "center", ha = "center")
    end
    line075black(); squareplot(); pretty_axes2()
    
end

