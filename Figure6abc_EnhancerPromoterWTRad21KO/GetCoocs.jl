suff = "1um_noz"

cooc = []

function get_coocs()

coocs = CSV.read(normpath("../CompleteSets/linkedlocus/", probe1*"__"*probe2*"__coocurrences"*suff*".csv"), DataFrames.DataFrame)
coocs[!,:BF_Enhancer] = coocs[!,"BF_"*probe2]
coocs[!,:BF_Gene] = coocs[!,"BF_"*probe1]
coocs[!,:Gene] = [probe1 for ii in 1:nrow(coocs)]
if probe1 == "Prdm1_intron"
        coocs[!,:Timepoint] = [split(ii, "_")[2] for ii in coocs[!,:Sample_Rep]]
        coocs = coocs[coocs[!,:Timepoint].!= "30", :]
        coocs = coocs[coocs[!,:Timepoint].!= "90", :]
        
    end

push!(cooc, coocs)
end

genefolder = "Il12bHSS1"
probe1= "Il12b_intron"
probe2 = "HSS1"

get_coocs()
genefolder = "Ifnb1L2"
probe1 = "Ifnb1forL2"
probe2 = "L2"

get_coocs()

genefolder = "Egr2"
probe1 = "Egr2_intron"
probe2 = "Egr2_enh"

get_coocs()

genefolder = "Prdm1"
probe1 = "Prdm1_intron"
probe2 = "Prdm1_enh"

get_coocs()

genefolder = "EhnPeliSertad"
probe2 = "Enh"
probe1 = "Peli1_intron"


get_coocs()

tb = join_in_all_common_columns(cooc...)
tb[!,:Genotype] = [split(ii, "_")[1] for ii in tb[!,:Sample_Rep]]
tb[!,:Timepoint] = [split(ii, "_")[2] for ii in tb[!,:Sample_Rep]]

tb[!,"log2_BF_Enhancer"] = log2.(tb[!,"BF_Enhancer"])
tb[!,"log2_BF_Gene"] = log2.(tb[!,"BF_Gene"]) 

CSV.write("EnhancerGEneBF.cvs", tb)