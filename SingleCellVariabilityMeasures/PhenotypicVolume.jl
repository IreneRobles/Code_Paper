
function scale(A)
    R"""a = scale($A)"""
    @rget a
end

function logPV(sceexp, cond; assay = "lnCPMplus1", permutations = 100, seed = 1223, cellssubsampled = 100,genesssubsampled = 1000)
    subsce = deepcopy(SingleCellExperiment.get_cells_with_this_characteristic(cond, :Sample, sceexp))
    dat = subsce.assays[assay]
    dat= transpose(dat) #feature on col, cell on row
    #remove non expressed genes
    bool = sum(dat, dims = 1).>0
    bool = reshape(bool,size(bool)[2])
    dat = dat[:,bool]
    dat = scale(dat)
    det = zeros(permutations)
    Random.seed!(seed)
    p = Progress(permutations, cond*": ")
    
   for ii in 1:permutations
        cellsample=shuffle(1:size(dat)[1])[1:cellssubsampled]
        genesample=shuffle(1:size(dat)[2])[1:genesssubsampled]
        dat = dat[cellsample,genesample]
        cdat = cov(dat)#calculate gene-gene covariance matrix using sampled genes
        svddat=svd(cdat) #singular value decomposition 
        eigen=svddat.S #get the eigen value 
        eigen=eigen[eigen.!=0]
        det[ii]= sum(log.(eigen)) #log sum of non-zero eigen values
        next!(p)
    end
    det
end

function PV_figure(table, col; u =5, h = maximum(table[!,col])+u, pvalue_ = :option2)
    sg = [split(ii, "_") for ii in table[!,:Sample]]
    table[!,:Genotype] = [ii[1] for ii in sg]
    table[!,:Timepoint] = [ii[2] for ii in sg]
    Seaborn.boxplot(data = Pandas.DataFrame(table), y = string(col), x = "Timepoint", hue = "Genotype", palette = ["darkgray", "red"], showfliers = false)
    
    tims = unique(table[!,:Timepoint])
    # Do test
    for tim in 1:length(tims)
        plt.plot([tim-1.25, tim-0.75], [h,h])
        subtb = table[table[!,:Timepoint].==tims[tim], :]
        wt = subtb[subtb[!,:Genotype].=="WT", col]
        rad = subtb[subtb[!,:Genotype].=="RAD21", col]
        
        p = 1 
        
        if pvalue_ == :option1
            p = minimum([(sum(wt .> rad)+1)/((nrow(subtb)/2)+1), (sum(wt .< rad)+1)/((nrow(subtb)/2)+1)])
        elseif pvalue_ == :option2
            p  = pvalue(MannWhitneyUTest(wt,rad))
        end
        println(p)
        x = (tim-1.25+tim-0.75)/2
        plt.annotate(round(p, sigdigits=2), xy = [x, h+u], va = "center", ha = "center")
    end
    line075black(); squareplot(); pretty_axes2()
    
end