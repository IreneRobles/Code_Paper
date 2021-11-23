
function timepoint(ii)
    if occursin("U 0", ii) 
        return "0"
    elseif occursin("U 1", ii) 
        return "0i"
    elseif occursin("U 2", ii) 
        return "8"
    elseif occursin("U 3", ii) 
        return "8i"
    else 
        "0"
    end
end


function genotype(ii)
    if occursin("V 0", ii) 
        return "WT"
    elseif occursin("V 1", ii) 
        return "Rad21KO"
    end
end


function genotype(ii)
    if occursin("V 0", ii) 
        return "WT"
    elseif occursin("V 1", ii) 
        return "Rad21KO"
    end
end


folder = "Fh1/"
folder = "../TSS_quantification/"*folder*"/"

gene = "Fh1"
e1 = get_data(folder, gene, 1, addrep = 1)
e2 = get_data(folder, gene, 2, addrep = 2)
e3 = get_data(folder, gene, 3, addrep = 3)
e4 = get_data(folder, gene, 4, addrep = 4)
exps = vcat(e1,e2, e3, e4)
exps[!,:Timepoint] = [timepoint(ii) for ii in exps[!,:Image]]
exps[!,:Genotype] = [genotype(ii) for ii in exps[!,:Image]]
exps[!,:Sample] = exps[!,:Genotype] .* "_" .* exps[!,:Timepoint]
CSV.write("GeneData/"*gene*".csv", exps)

folder = "Hprt/"
folder = "../TSS_quantification/"*folder*"/"

gene = "Hprt"


e1 = get_data(folder, gene, 1, addrep = 1)
e2 = get_data(folder, gene, 2, addrep = 2)
e3 = get_data(folder, gene, 3, addrep = 3)
e4 = get_data(folder, gene, 4, addrep = 4)
exps = vcat(e1,e2, e3, e4)

exps[!,:Timepoint] = [timepoint(ii) for ii in exps[!,:Image]]
exps[!,:Genotype] = [genotype(ii) for ii in exps[!,:Image]]
exps[!,:Sample] = exps[!,:Genotype] .* "_" .* exps[!,:Timepoint]
CSV.write("GeneData/"*gene*".csv", exps)


