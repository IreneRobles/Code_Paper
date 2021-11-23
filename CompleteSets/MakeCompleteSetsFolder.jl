colstosave = [ "Image","Cell", "Genotype", "Timepoint", "Rep","TSS1_r2", "TSS2_r2"]


gene = "Il12b_intron_original"
set1 = CSV.read("GeneData/"*gene*".csv", DataFrame)

gene = "Ifit1"
ifit1 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Ifit1.csv", ifit1[!, colstosave])




colstosave = [ "Image","Cell", "Genotype", "Timepoint", "Rep","TSS1_r2", "TSS2_r2"]
gene = "Ifnb1_original"
set1 = CSV.read("GeneData/"*gene*".csv", DataFrame)
gene = "Ifnb1forL2"
set2 = CSV.read("GeneData/"*gene*".csv", DataFrame)

add_to_rep = maximum(unique(set1[!,:Rep]))
set2[!, :Rep] = set2[!, :Rep] .+ add_to_rep
ifnb1 = join_in_all_common_columns(set1, set2)
ifnb1 = ifnb1[ifnb1[!,:TSS1_r2].!="NA", :]
ifnb1 = ifnb1[ifnb1[!,:TSS2_r2].!="NA", :]
CSV.write("CompleteSets/Ifnb1.csv", ifnb1[!, colstosave])


colstosave2 = [ "Image","Cell", "Genotype", "Timepoint", "Rep","TSS1_r2", "TSS2_r2", "N_exon"]
CSV.write("CompleteSets/Ifnb1_mature.csv", set1[!, colstosave2])



gene = "Il12b_intron_original"
set1 = CSV.read("GeneData/"*gene*".csv", DataFrame)

colstosave2 = [ "Image","Cell", "Genotype", "Timepoint", "Rep","TSS1_r2", "TSS2_r2", "N_exon"]
CSV.write("CompleteSets/Il12b_mature_nascent.csv", set1[!, colstosave2])

gene = "Il12b_intron"
set2 = CSV.read("GeneData/"*gene*".csv", DataFrame)
add_to_rep = maximum(unique(set1[!,:Rep]))
set2[!, :Rep] = set2[!, :Rep] .+ add_to_rep
il12b = join_in_all_common_columns(set1, set2)
il12b = il12b[il12b[!,:TSS1_r2].!="NA", :]
il12b = il12b[il12b[!,:TSS2_r2].!="NA", :]
CSV.write("CompleteSets/Il12b.csv", il12b[!, colstosave])



gene = "Egr2_intron"
egr2 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Egr2.csv", egr2[!, colstosave])


gene = "Prdm1_intron"
prdm1 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Prdm1.csv", prdm1[!, colstosave])


gene = "Fh1"
fh1 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Fh1.csv", fh1[!, colstosave])


gene = "Hprt"
hprt = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Hprt.csv", hprt[!, colstosave])

gene = "Ifit1"
ifit1 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Ifit1.csv", ifit1[!, colstosave])
colstosave2 = [ "Image","Cell", "Genotype", "Timepoint", "Rep","TSS1_r2", "TSS2_r2", "N_exon"]
CSV.write("CompleteSets/Ifit1_mature.csv", set1[!, colstosave2])

gene = "Cxcl10"
cxcl10 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Cxcl10.csv", cxcl10[!, colstosave])

gene = "Cxcl10"
Cxcl10 = CSV.read("GeneData/"*gene*".csv", DataFrames.DataFrame)
Cxcl101 = Cxcl10[Cxcl10[!,:Genotype].=="WT", :]
Cxcl102 = Cxcl10[Cxcl10[!,:Genotype].=="Rad21KO", :]
Cxcl10 = vcat(Cxcl101, Cxcl102)
Cxcl10 = Cxcl10[Cxcl10[!,:TSS1_r2].!="NA", :]
Cxcl10 = Cxcl10[Cxcl10[!,:TSS2_r2].!="NA", :]
rename!(Cxcl10,"N_thres_Total" =>"N_exon",
    "N_thres_Nuc" =>"N_exon_Nuc",)
CSV.write("CompleteSets/Cxcl10_mature.csv", Cxcl10[!, colstosave2])
Cxcl10 = CSV.read("CompleteSets/Cxcl10_mature.csv", DataFrames.DataFrame)

gene = "Peli1_intron_original"
set1 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Peli1_mature_nascent.csv", set1[!, colstosave2])
gene = "Peli1_intron"
set2 = CSV.read("GeneData/"*gene*".csv", DataFrame)

add_to_rep = maximum(unique(set1[!,:Rep]))
set2[!, :Rep] = set2[!, :Rep] .+ add_to_rep
peli1 = join_in_all_common_columns(set1, set2)
CSV.write("CompleteSets/Peli1.csv", peli1[!, colstosave])

gene = "Sertad2"
set2 = CSV.read("GeneData/Sertad2_intron.csv", DataFrame)
sertad = set2
CSV.write("CompleteSets/Sertad2.csv", sertad[!, colstosave])



gene = "Egr2_enh"
egr2 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Egr2_enh.csv", egr2[!, colstosave])


gene = "Prdm1_enh"
prdm1 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Prdm1_enh.csv", prdm1[!, colstosave])


gene = "Enh"
fh1 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/Enh.csv", fh1[!, colstosave])


gene = "HSS1"
hprt = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/HSS1.csv", hprt[!, colstosave])

gene = "L2"
ifit1 = CSV.read("GeneData/"*gene*".csv", DataFrame)
CSV.write("CompleteSets/L2.csv", ifit1[!, colstosave])