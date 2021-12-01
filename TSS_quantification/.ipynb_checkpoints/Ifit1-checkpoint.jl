if !in("Ifit1", readdir())
    mkdir("Ifit1")
end

ENV["Code"] = "../../Code"
for folder in readdir(ENV["Code"]); push!(LOAD_PATH, normpath(ENV["Code"], folder));end

include("../../Code/TSSs/src/TSSs.jl")
using ProgressMeter
using DataFrames
using Distances
using CSV
#using RCall
using NoLongerProblems_FileHandling
using NoLongerProblems
using NoLongerProblems_Pandas
using DataFrames
using FQfiles
using HypothesisTests
using MultipleTesting
import Pandas
using PyPlot
using PrettyPlotting


function CellInfo(dir, exp_df)
    image_file = get_files_ending_with(dir, "Image.csv")
     cell_file = get_files_ending_with(dir, "TrueCells.csv")

    if .&(length(cell_file) == 1,  length(image_file) == 1) # Make sure that there is only one file with the image information
        image_file = normpath(dir, image_file[1])
        cell_file = normpath(dir, cell_file[1])
    else
        error("Image.csv or Cell file not found")
    end
    
    ima = DataFrames.DataFrame(CSV.read(image_file, DataFrames.DataFrame))[!, [:FileName_DAPI, :ImageNumber]]
    ima[!,:Image] = [split(split(ii, "_C1")[1], "_MAX")[1] for ii in ima[!,:FileName_DAPI]]
    
    ima[!,:WELL] = [split(ii, " ") for ii in ima[!,:Image]]
    ima[!,:WELL] = [1 for ii in ima[!,:WELL]]
    
     ima[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in ima[!,:FileName_DAPI]]
    
    exp_df[!,:Well] = [split(ii, " (")[1] for ii in exp_df[!,:Well]]

    im_ = innerjoin(exp_df, ima, on = :Well, makeunique=true)
        

    
   
     cells = DataFrames.DataFrame(CSV.read(cell_file, DataFrames.DataFrame))
    cells = innerjoin(im_, cells, on = :ImageNumber)
    
    cells[!,:Image_Cell] = [cells[ii, :Image]*"__Cell_CP_"*string(cells[ii, :ObjectNumber]) for ii in 1:nrow(cells)]
    return dropmissing(cells)
    
end



exp = CSV.read("IbidiChambers/IbidiChamberSlide_Exp1.csv",  DataFrames.DataFrame)
exp[!,:ChamberSlide] = [if split(ii, "(")[2][1] == '1' "UT" else "3H" end for ii in exp[!,:Well]]
exp[!,:WELL] = ["1" for ii in exp[!,:Well]]
exp
cp_dir = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp1/CP_results"
cells1 = CellInfo(cp_dir, add_probetypes(exp));
cp_dir = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp2/quick_analysis_qq/CP_results"
cells2 = CellInfo(cp_dir, add_probetypes(exp));
cp_dir = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp3/CP_results"
cells3 = CellInfo(cp_dir, add_probetypes(exp));


# Define where are the outlines

root1 = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp1"
root2 = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp2/quick_analysis_qq"
root3 = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp3/"


tss_c4 = root1 * "/Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 4; xy = 0.189, zx = 0.5)
if !in("TSS_raw", readdir())
    mkdir("TSS_raw")
end
CSV.write("TSS_raw/Ifit1_exp1.csv", d)



tss_c4 = root2 * "/Segmentation_C4/_FQ_outline/_TS_detect"
imagesfolder = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp2/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d= TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 4; xy = 0.189, zx = 0.5)
if !in("TSS_raw", readdir())
    mkdir("TSS_raw")
end
CSV.write("TSS_raw/Ifit1_exp2.csv", d)

# Define where are the outlines


tss_c4 = root3 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root3 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d= TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 4; xy = 0.189, zx = 0.5)
if !in("TSS_raw", readdir())
    mkdir("TSS_raw")
end
CSV.write("TSS_raw/Ifit1_exp3.csv", d)

if !in("TSS_avgdot", readdir())
    mkdir("TSS_avgdot")
end
root2 = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp2/quick_analysis_qq"

tss2 = CSV.read("TSS_raw/Ifit1_exp1.csv", DataFrames.DataFrame)
tss3 = CSV.read("TSS_raw/Ifit1_exp2.csv", DataFrames.DataFrame)
tss4 = CSV.read("TSS_raw/Ifit1_exp3.csv", DataFrames.DataFrame)


dot2_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type1_TSS/_mRNA_AVG_ns.tif"); radious = 1)*nor
dot3_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type1_TSS/_mRNA_AVG_ns.tif"); radious = 1)*nor
dot4_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type1_TSS/_mRNA_AVG_ns.tif"); radious = 1)*nor

dot2_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type1_TSS/_mRNA_AVG_ns.tif"); radious = 2)*nor
dot3_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type1_TSS/_mRNA_AVG_ns.tif"); radious = 2)*nor
dot4_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type1_TSS/_mRNA_AVG_ns.tif"); radious = 2)*nor

dot2_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type1_TSS/_mRNA_AVG_ns.tif"); radious = 3)*nor
dot3_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type1_TSS/_mRNA_AVG_ns.tif"); radious = 3)*nor
dot4_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type1_TSS/_mRNA_AVG_ns.tif"); radious = 3)*nor


tss2[!,:TSS1_r1] = tss2[!,:locus1_int1_TSS2] ./ dot2_r1
tss2[!,:TSS1_r2] = tss2[!,:locus1_int2_TSS2] ./ dot2_r2
tss2[!,:TSS1_r3] = tss2[!,:locus1_int3_TSS2] ./ dot2_r3


tss3[!,:TSS1_r1] = tss3[!,:locus1_int1_TSS2] ./ dot3_r1
tss3[!,:TSS1_r2] = tss3[!,:locus1_int2_TSS2] ./ dot3_r2
tss3[!,:TSS1_r3] = tss3[!,:locus1_int3_TSS2] ./ dot3_r3


tss4[!,:TSS1_r1] = tss4[!,:locus1_int1_TSS2] ./ dot4_r1
tss4[!,:TSS1_r2] = tss4[!,:locus1_int2_TSS2] ./ dot4_r2
tss4[!,:TSS1_r3] = tss4[!,:locus1_int3_TSS2] ./ dot4_r3


tss2[!,:TSS2_r1] = tss2[!,:locus2_int1_TSS2] ./ dot2_r1
tss2[!,:TSS2_r2] = tss2[!,:locus2_int2_TSS2] ./ dot2_r2
tss2[!,:TSS2_r3] = tss2[!,:locus2_int3_TSS2] ./ dot2_r3


tss3[!,:TSS2_r1] = tss3[!,:locus2_int1_TSS2] ./ dot3_r1
tss3[!,:TSS2_r2] = tss3[!,:locus2_int2_TSS2] ./ dot3_r2
tss3[!,:TSS2_r3] = tss3[!,:locus2_int3_TSS2] ./ dot3_r3


tss4[!,:TSS2_r1] = tss4[!,:locus2_int1_TSS2] ./ dot4_r1
tss4[!,:TSS2_r2] = tss4[!,:locus2_int2_TSS2] ./ dot4_r2
tss4[!,:TSS2_r3] = tss4[!,:locus2_int3_TSS2] ./ dot4_r3


CSV.write("TSS_avgdot/Ifit1_exp1.csv", tss2)
CSV.write("TSS_avgdot/Ifit1_exp2.csv", tss3)
CSV.write("TSS_avgdot/Ifit1_exp3.csv", tss4)

exp = CSV.read("IbidiChambers/IbidiChamberSlide_Ifit1_mixexps.csv", DataFrames.DataFrame)
exp[!,:ChamberSlide] = [if split(ii, "(")[2][1] == '1' "UT" else "3H" end for ii in exp[!,:Well]]
exp[!,:WELL] = ["1" for ii in exp[!,:Well]]
exp[!,:Genotype] = exp[!,:Cells]
exp[!,:Well] = [split(ii, " (")[1] for ii in exp[!,:Well]]


cp_dir = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp1/CP_results"
cells1 = CellInfo(cp_dir, add_probetypes(exp));
cp_dir = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp2/quick_analysis_qq/CP_results"
cells2 = CellInfo(cp_dir, add_probetypes(exp));
cp_dir = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp3/CP_results"
cells3 = CellInfo(cp_dir, add_probetypes(exp));
exp1 = CSV.read("TSS_avgdot/Ifit1_exp1.csv", DataFrames.DataFrame)
exp1[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in exp1[!,:Image]]
exp1 = innerjoin(exp1, exp, on= "Well")
CSV.write("Ifit1/Ifit1_exp1_TSS.csv", exp1)
CSV.write("Ifit1/Ifit1_exp1_CP.csv", cells1)
exp1 = CSV.read("TSS_avgdot/Ifit1_exp2.csv", DataFrames.DataFrame)
exp1[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in exp1[!,:Image]]
exp1 = innerjoin(exp1, exp, on= "Well")
CSV.write("Ifit1/Ifit1_exp2_TSS.csv", exp1)
CSV.write("Ifit1/Ifit1_exp2_CP.csv", cells2)
exp1 = CSV.read("TSS_avgdot/Ifit1_exp3.csv", DataFrames.DataFrame)
exp1[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in exp1[!,:Image]]
exp1 = innerjoin(exp1, exp, on= "Well")
CSV.write("Ifit1/Ifit1_exp3_TSS.csv", exp1)
CSV.write("Ifit1/Ifit1_exp3_CP.csv", cells3)

fq_dir1 = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp1/type1_TSS/"
fq_dir2 = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp2/quick_analysis_qq/type1_TSS/"
fq_dir3 = "/Volumes/lymphdev\$/SarahWells/Confocal/CellMix_Rad21KO_WT_Ifit1_Mx1/Exp3/type1_TSS/"


fq1 = FQ_summary_MATURE(fq_dir1)
rename!(fq1, :N_thres_Total => :N_exon)
rename!(fq1, :N_thres_Nuc => :N_exon_Nuc)
fq2 = FQ_summary_MATURE(fq_dir2)
rename!(fq2, :N_thres_Total => :N_exon)
rename!(fq2, :N_thres_Nuc => :N_exon_Nuc)
fq3 = FQ_summary_MATURE(fq_dir3)
rename!(fq3, :N_thres_Total => :N_exon)
rename!(fq3, :N_thres_Nuc => :N_exon_Nuc)

CSV.write("Ifit1/Ifit1_exp1_FQ.csv", fq1)
CSV.write("Ifit1/Ifit1_exp2_FQ.csv", fq2)
CSV.write("Ifit1/Ifit1_exp3_FQ.csv", fq3)
