BF = Dict()

bf = calculate_bf_by_rep(get_completeset("Peli1"), "TSS", "Peli1",limit = 1)
# Only consider 1 timepoint
bf[!,:Timepoint] = [if ii .== "135" "120" else ii end for ii in bf[!,:Timepoint]]
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "120"
bool = bool1.| bool2
BF["Peli1"] = bf[bool, :]

bf = calculate_bf_by_rep(get_completeset("Il12b"), "TSS", "Il12b",limit = 1)
# Only consider 1 timepoint
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "90"
bool = bool1.| bool2
BF["Il12b"] = bf[bool, :]


bf = calculate_bf_by_rep(get_completeset("Ifnb1"), "TSS", "Ifnb1", limit = 1)
# Only consider 1 timepoint
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "90"

bool3 = bf[!,:Rep] .!= 3

bool = (bool1.| bool2).*bool3
BF["Ifnb1"] = bf[bool, :]


bf = calculate_bf_by_rep(get_completeset("Ifit1"), "TSS", "Ifit1",limit = 1)
# Only consider 1 timepoint
BF["Ifit1"] = bf

bf = calculate_bf_by_rep(get_completeset("Cxcl10"), "TSS", "Cxcl10",limit = 1)
# Only consider 1 timepoint
BF["Cxcl10"] = bf

bf = calculate_bf_by_rep(get_completeset("Sertad2"), "TSS", "Sertad2",limit = 1)
# Only consider 1 timepoint
BF["Sertad2"] = bf


bf = calculate_bf_by_rep(get_completeset("Prdm1"), "TSS", "Prdm1",limit = 0)
# Only consider 1 timepoint
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "60"
bool = bool1.| bool2
BF["Prdm1"] = bf[bool, :]

bf = calculate_bf_by_rep(get_completeset("Egr2"), "TSS", "Egr2",limit = 1)
# Only consider 1 timepoint
BF["Egr2"] = bf


bf = calculate_bf_by_rep(get_completeset("Fh1"), "TSS", "Fh1",limit = 1)
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "8"
bool = bool1.| bool2
# Only consider 1 timepoint
BF["Fh1"] = bf[bool, :]


bf = calculate_bf_by_rep(get_completeset("Hprt"), "TSS", "Hprt",limit = 1)
bool1 = bf[!,:Timepoint] .== "0" 
bool2 = bf[!,:Timepoint] .== "8"
bool = bool1.| bool2
# Only consider 1 timepoint
BF["Hprt"] = bf[bool, :]

bf = calculate_bf_by_rep(get_genedata("Enh"), "TSS", "Enh", limit = 1)
BF["Enh"] = bf
bf = calculate_bf_by_rep(get_genedata("L2"), "TSS", "L2", limit = 1)
BF["L2"] = bf
bf = calculate_bf_by_rep(get_genedata("Egr2_enh"), "TSS", "Egr2_enh", limit = 1)
BF["Egr2_enh"] = bf
bf = calculate_bf_by_rep(get_genedata("Prdm1_enh"), "TSS", "Prdm1_enh", limit =1)
bf = bf[bf[!,:Timepoint].!= "30", :]
bf = bf[bf[!,:Timepoint].!= "90", :]
BF["Prdm1_enh"] = bf
bf = calculate_bf_by_rep(get_genedata("HSS1"), "TSS", "HSS1", limit = 0)
BF["HSS1"] = bf
BFs = join_in_all_common_columns([BF[key] for key in keys(BF)]...)

CSV.write("BurstFrequencies_genes_enhancers.csv", BFs)