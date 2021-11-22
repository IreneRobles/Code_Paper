# Il12b
gene = "Il12b_intron_original"
folder = "../TSS_quantification/Il12b_original/"
exps = []
for expnumber in ["1", "2", "3", "4"]
exp = CSV.read(folder* gene*"_exp"*expnumber*"_FQ.csv", DataFrames.DataFrame)
exp = column_fusion(exp, :Image, :Cell)
t = CSV.read(folder* gene*"_exp"*expnumber*"_TSS.csv", DataFrames.DataFrame)
t[!,:Image_Cell] = t[!,:Image] .*"__".*string.(t[!,:Cell])
exp = innerjoin(exp, t, on = :Image_Cell, makeunique=true)
rename!(exp, :N_thres_Total => :N_intron, :N_thres_Nuc => :N_intron_Nuc;)
mature = CSV.read(folder* gene*"_exp"*expnumber*"_FQ_mature.csv", DataFrames.DataFrame)
rename!(mature, :N_thres_Total => :N_exon, :N_thres_Nuc => :N_exon_Nuc;)
exp = innerjoin(exp, mature[!,["N_exon","N_exon_Nuc","Image_Cell"]], on = :Image_Cell)
exp[!, :Rep] = [expnumber for ii in 1:nrow(exp)]
    cp = CSV.read(folder * gene*"_exp"*expnumber*"_CP.csv", DataFrames.DataFrame)[!,["Genotype","Timepoint","Image_Cell"]]
    cp[!,"Image_Cell"] = [replace(ii, "MAX_" => "") for ii in cp[!,"Image_Cell"]]
    exp = sample_from_genotype_timepoint(sort!(innerjoin(exp, cp, on = :Image_Cell, makeunique=true), :Genotype, rev = true))
push!(exps, exp)
end
exps = vcat(exps...)
exps[!,:TSS1_r2] = [parse(Float64, split(split(ii, "(")[end], ")")[1]) for ii in  exps[!,:TSS1_r2]]
exps[!,:TSS2_r2] = [parse(Float64, split(split(ii, "(")[end], ")")[1]) for ii in  exps[!,:TSS2_r2]]
CSV.write("GeneData/Il12b_intron_original.csv", exps)

# Peli1
gene = "Peli1_intron_original"
folder = "../TSS_quantification/Peli1_original/"
exps = []
for expnumber in ["1", "2", "3"]
exp = CSV.read(folder* gene*"_exp"*expnumber*"_FQ.csv", DataFrames.DataFrame)
exp = column_fusion(exp, :Image, :Cell)
t = CSV.read(folder* gene*"_exp"*expnumber*"_TSS.csv", DataFrames.DataFrame)
t[!,:Image_Cell] = t[!,:Image] .*"__".*string.(t[!,:Cell])
exp = innerjoin(exp, t, on = :Image_Cell, makeunique=true)
rename!(exp, :N_thres_Total => :N_intron, :N_thres_Nuc => :N_intron_Nuc;)
mature = CSV.read(folder* gene*"_exp"*expnumber*"_FQ_mature.csv", DataFrames.DataFrame)
rename!(mature, :N_thres_Total => :N_exon, :N_thres_Nuc => :N_exon_Nuc;)
exp = innerjoin(exp, mature[!,["N_exon","N_exon_Nuc","Image_Cell"]], on = :Image_Cell)
exp[!, :Rep] = [expnumber for ii in 1:nrow(exp)]
    cp = CSV.read(folder * gene*"_exp"*expnumber*"_CP.csv", DataFrames.DataFrame)[!,["Genotype","Timepoint","Image_Cell"]]
    cp[!,"Image_Cell"] = [replace(ii, "MAX_" => "") for ii in cp[!,"Image_Cell"]]
    exp = sample_from_genotype_timepoint(sort!(innerjoin(exp, cp, on = :Image_Cell, makeunique=true), :Genotype, rev = true))
push!(exps, exp)
end
exps = vcat(exps...)
exps[!,:TSS1_r2] = [parse(Float64, split(split(ii, "(")[end], ")")[1]) for ii in  exps[!,:TSS1_r2]]
exps[!,:TSS2_r2] = [parse(Float64, split(split(ii, "(")[end], ")")[1]) for ii in  exps[!,:TSS2_r2]]
CSV.write("GeneData/Peli1_intron_original.csv", exps)

# Ifnb1
gene = "Ifnb1_original"
folder = "../TSS_quantification/Ifnb1_original/"
exps = []
for expnumber in ["1", "2",]
exp = CSV.read(folder* gene*"_exp"*expnumber*"_FQ.csv", DataFrames.DataFrame)
exp = column_fusion(exp, :Image, :Cell)
t = CSV.read(folder* gene*"_exp"*expnumber*"_TSS.csv", DataFrames.DataFrame)
t[!,:Image_Cell] = t[!,:Image] .*"__".*string.(t[!,:Cell])
exp = innerjoin(exp, t, on = :Image_Cell, makeunique=true)
rename!(exp, :N_thres_Total => :N_exon, :N_thres_Nuc => :N_exon_Nuc;)
cp = CSV.read(folder * gene*"_exp"*expnumber*"_CP.csv", DataFrames.DataFrame)[!,["Genotype","Timepoint","Image_Cell"]]
    cp[!,"Image_Cell"] = [replace(ii, "MAX_" => "") for ii in cp[!,"Image_Cell"]]
    exp = sample_from_genotype_timepoint(sort!(innerjoin(exp, cp, on = :Image_Cell, makeunique=true), :Genotype, rev = true))
push!(exps, exp)
end
exps = vcat(exps...)
exps[!,:TSS1_r2] = [parse(Float64, split(split(ii, "(")[end], ")")[1]) for ii in  exps[!,:TSS1_r2]]
exps[!,:TSS2_r2] = [parse(Float64, split(split(ii, "(")[end], ")")[1]) for ii in  exps[!,:TSS2_r2]]

CSV.write("GeneData/Ifnb1_original.csv", exps)