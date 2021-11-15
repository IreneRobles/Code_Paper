if !in("Il12bHSS1", readdir())
    mkdir("Il12bHSS1")
end

function normalise_avgdot(df, genename, expn, typeprobe, ibidislide; root = pwd(), genefolder = "Gene", cellinfo = CellInfo)
    dot1_r1 = TSSs.int_brightest_pixel(TSSs.read_tiff_as_gray(root * typeprobe *"/_mRNA_AVG_ns.tif").*nor; radious = 1)
    dot1_r2 = TSSs.int_brightest_pixel(TSSs.read_tiff_as_gray(root * typeprobe * "/_mRNA_AVG_ns.tif").*nor; radious = 2)
    dot1_r3 = TSSs.int_brightest_pixel(TSSs.read_tiff_as_gray(root * typeprobe * "/_mRNA_AVG_ns.tif").*nor; radious = 3)

    df[!,:TSS1_r1] = df[!,:locus1_int1_TSS2] ./ dot1_r1
    df[!,:TSS1_r2] = df[!,:locus1_int2_TSS2] ./ dot1_r2
    df[!,:TSS1_r3] = df[!,:locus1_int3_TSS2] ./ dot1_r3
    df[!,:TSS2_r1] = df[!,:locus2_int1_TSS2] ./ dot1_r1
    df[!,:TSS2_r2] = df[!,:locus2_int2_TSS2] ./ dot1_r2
    df[!,:TSS2_r3] = df[!,:locus2_int3_TSS2] ./ dot1_r3
    
    
    CSV.write("TSS_avgdot/"*genename*"_exp"*"$expn"*".csv", df)
    ibidislide[!,:Well] = [split(ii, " (")[1] for ii in ibidislide[!,:Well]]
    cp_dir = root1*"CP_results"
     cells1 = cellinfo(cp_dir, add_probetypes(exp));
    exp1 = CSV.read("TSS_avgdot/"*gene*"_exp"*"$expn"*".csv", DataFrame)
      well = []
    for ii in exp1[!,:Image]
        try
            push!(well,split(split(ii, "S 0_")[2], "_X")[1])
        catch
            push!(well,split(split(ii, "_Pos")[1], "1_0_")[2])
        end
    end
    exp1[!,:Well] = well
    exp1 = innerjoin(exp1, ibidislide, on= "Well")
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_TSS.csv", exp1)
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_CP.csv", cells1)
    fq_dir1 = root1*typeprobe*"/"
    fq1 = FQ_summary_MATURE(fq_dir1)
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ.csv", fq1)
end

function CellInfo1(dir, exp_df)
    image_file = get_files_ending_with(dir, "Image.csv")
     cell_file = get_files_ending_with(dir, "TrueCells.csv")

    if .&(length(cell_file) == 1,  length(image_file) == 1) # Make sure that there is only one file with the image information
        image_file = normpath(dir, image_file[1])
        cell_file = normpath(dir, cell_file[1])
    else
        error("Image.csv or Cell file not found")
    end
    
    ima = DataFrames.DataFrame(CSV.read(image_file, DataFrame))[:, [:FileName_DAPI, :ImageNumber]]
    ima[!,:Image] = [split(split(ii, "_C1")[1], "_MAX")[1] for ii in ima[!,:FileName_DAPI]]
    
    ima[!,:WELL] = [split(ii, " ") for ii in ima[!,:Image]]
    ima[!,:WELL] = [1 for ii in ima[!,:WELL]]
    well = []
    for ii in ima[!,:FileName_DAPI]
        try
            push!(well,split(split(ii, "S 0_")[2], "_X")[1])
        catch
            push!(well,split(split(ii, "_Pos")[1], "1_0_")[2])
        end
    end
    
     ima[!,:Well] = well
    
    exp_df[!,:Well] = [split(ii, " (")[1] for ii in exp_df[!,:Well]]

    im_ = innerjoin(exp_df, ima, on = :Well, makeunique=true)
    cells = DataFrames.DataFrame(CSV.read(cell_file, DataFrame))
    cells = outerjoin(im_, cells, on = :ImageNumber)
    
    cells[!,:Image_Cell] = [cells[ii, :Image]*"__Cell_CP_"*string(cells[ii, :ObjectNumber]) for ii in 1:nrow(cells)]
    return dropmissing(cells)
    
end

function removeimages(df, images...)
    splitname = split.(df[!,:Image], "_")
    imnum = [parse(Int, ii[end]) for ii in splitname]
    imbool = .! in.(imnum,images)
    df[imbool, :]
end



ROOT1 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Il12bBD/Rep1/"
ROOT2 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Il12bBD/Rep2/"
ROOT3 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Il12bBD/Rep5/"
ROOT4 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Il12bBD/Rep6/"
ROOT5 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Il12bBD/Rep78/"
ROOT6 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Il12bBD/Rep9/"
ROOT7 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Il12bBD/Rep10/"

root1 = ROOT1
expn = "1"

gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT2
expn = "2"

gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

root1 = ROOT3
expn = "3"

gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

root1 = ROOT4
expn = "4"

gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT4
expn = "4"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect40"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

root1 = ROOT5
expn = "5_30"

gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect30"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)


root1 = ROOT5
expn = "5_35"

gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect35"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)


root1 = ROOT5
expn = "5_45"

gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect45"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT5
expn = "5_40"

gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect40"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT5
expn = "5_37"

gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

root1 = ROOT6
expn = "6_60"

gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect_60"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

expn = "6_70"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect_70"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

expn = "6_100"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect_100"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

root1 = ROOT6

expn = "6_60"
gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect_60"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

expn = "6_65"
gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect_65"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

expn = "6_70"
gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect_70"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

expn = "6_75"
gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect_75"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

expn = "6_50"
gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect_50"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT7

expn = "7_40"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect_40"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

expn = "7_50"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect_50"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

expn = "7_60"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect_60"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)


expn = "7_35"
gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect_35"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

expn = "7_47"
gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect_47"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

expn = "7_60"
gene = "HSS1_BD"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect_60"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)


folder = "Il12bHSS1BD"
root1 = ROOT1
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp1.csv", DataFrame)
expn = 1
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)
folder = "Il12bHSS1BD"
root1 = ROOT2
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp2.csv", DataFrame)
expn = 2
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")

folder = "Il12bHSS1BD"
root1 = ROOT3
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp3.csv", DataFrame)
expn = 3
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")




folder = "Il12bHSS1BD"
root1 = ROOT4
ims = [136,312,166,310,240,308,242,296,520,440,298, 484]
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp4.csv", DataFrame)
expn = 4
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
tss1 = removeimages(tss1, ims)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
tss1 = removeimages(tss1, ims)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")


folder = "Il12bHSS1BD"
root1 = ROOT1
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp1.csv", DataFrame)
expn = 1
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)
folder = "Il12bHSS1BD"
root1 = ROOT2
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp2.csv", DataFrame)
expn = 2
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")

folder = "Il12bHSS1BD"
root1 = ROOT3
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp3.csv", DataFrame)
expn = 3
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")




folder = "Il12bHSS1BD"
root1 = ROOT4
ims = [136,312,166,310,240,308,242,296,520,440,298, 484]
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp4.csv", DataFrame)
expn = 4
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
tss1 = removeimages(tss1, ims)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
tss1 = removeimages(tss1, ims)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")

function normalise_avgdot2(df, genename, expn, typeprobe, ibidislide; root = pwd(), genefolder = "Gene", cellinfo = CellInfo)
    dot1_r1 = TSSs.int_brightest_pixel(TSSs.read_tiff_as_gray(root * typeprobe *"/_mRNA_AVG_ns.tif").*nor; radious = 1)
    dot1_r2 = TSSs.int_brightest_pixel(TSSs.read_tiff_as_gray(root * typeprobe * "/_mRNA_AVG_ns.tif").*nor; radious = 2)
    dot1_r3 = TSSs.int_brightest_pixel(TSSs.read_tiff_as_gray(root * typeprobe * "/_mRNA_AVG_ns.tif").*nor; radious = 3)

    df[!,:TSS1_r1] = df[!,:locus1_int1_TSS2] ./ dot1_r1
    df[!,:TSS1_r2] = df[!,:locus1_int2_TSS2] ./ dot1_r2
    df[!,:TSS1_r3] = df[!,:locus1_int3_TSS2] ./ dot1_r3
    df[!,:TSS2_r1] = df[!,:locus2_int1_TSS2] ./ dot1_r1
    df[!,:TSS2_r2] = df[!,:locus2_int2_TSS2] ./ dot1_r2
    df[!,:TSS2_r3] = df[!,:locus2_int3_TSS2] ./ dot1_r3
    
    
    CSV.write("TSS_avgdot/"*genename*"_exp"*"$expn"*".csv", df)
    ibidislide[!,:Well] = [split(ii, " (")[1] for ii in ibidislide[!,:Well]]
    cp_dir = root1*"CP_results"
     cells1 = cellinfo(cp_dir, add_probetypes(exp));
    exp1 = CSV.read("TSS_avgdot/"*gene*"_exp"*"$expn"*".csv", DataFrame)
    exp1[!,:Slide] = [ii[1:6] for ii in exp1[!,:Image]]
      well = []
    for ii in exp1[!,:Image]
        try
            push!(well,split(split(ii, "S 0_")[2], "_X")[1])
        catch
            push!(well,split(split(ii, "_Pos")[1], "1_0_")[2])
        end
    end
    exp1[!,:Well] = well
    exp1 = innerjoin(exp1, ibidislide, on= ["Well", "Slide"])
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_TSS.csv", exp1)
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_CP.csv", cells1)
    fq_dir1 = root1*typeprobe*"/"
    fq1 = FQ_summary_MATURE(fq_dir1)
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ.csv", fq1)
end


folder = "Il12bHSS1BD"
root1 = ROOT6
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp9.csv", DataFrame)
expn = "6_60"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "6_70"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "6_100"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "6_50"
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "6_60"
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "6_65"
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "6_70"
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "6_75"
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

folder = "Il12bHSS1BD"
root1 = ROOT7
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp9.csv", DataFrame)
expn = "7_40"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "7_50"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "7_60"
gene = "Il12b_intron_BD"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "7_35"
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "7_47"
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

expn = "7_60"
gene = "HSS1_BD"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1",  cellinfo = CellInfo1)

folder = "Il12bHSS1BD"
root1 = ROOT5
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_BD_exp78.csv", DataFrame)
gene = "HSS1_BD"
typeprobe = "type1"
expn = "5_30"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot2(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
expn = "5_35"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot2(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
expn = "5_37"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot2(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
expn = "5_40"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot2(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
expn = "5_45"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot2(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")


