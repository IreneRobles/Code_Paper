function normalise_avgdot(df, genename, expn, typeprobe, ibidislide; root = pwd(), genefolder = "Gene")
    
    if !in(genefolder, readdir())
        mkdir(genefolder)
    end
    
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
    cells1 = CellInfo(cp_dir, add_probetypes(exp));
    exp1 = CSV.read("TSS_avgdot/"*gene*"_exp"*"$expn"*".csv", DataFrame)
    exp1[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in exp1[!,:Image]]
    exp1 = innerjoin(exp1, ibidislide, on= "Well")
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_TSS.csv", exp1)
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_CP.csv", cells1)
    fq_dir1 = root1*typeprobe*"/"
    fq1 = FQ_summary_MATURE(fq_dir1)
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ.csv", fq1)
end

root1 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Egr2/Exp1/"
expn = 1
gene = "Egr2_enh"
typeprobe = "type6"
tss_c4 = root1 * "Segmentation_type6/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

gene = "Egr2_intron"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

root1 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Egr2/Exp2/"
expn = 2
gene = "Egr2_enh"
typeprobe = "type6"
tss_c4 = root1 * "Segmentation_type6/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

gene = "Egr2_intron"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)


root1 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Egr2/Exp1/"
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Egr2_20210125.csv", DataFrame)
expn = 1
gene = "Egr2_enh"
typeprobe = "type6"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Egr2")


gene = "Egr2_intron"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Egr2")

root1 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Egr2/Exp2/"
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Egr2_20210205.csv", DataFrame)
expn = 2
gene = "Egr2_enh"
typeprobe = "type6"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Egr2")

gene = "Egr2_intron"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Egr2")

ROOT1 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Ifnb1/Exp1/"

root1 = ROOT1

expn = 1
gene = "Ifnb1forL2"
typeprobe = "type4"
channel= 3
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)

d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, channel; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

gene = "L2"
typeprobe = "type6"
channel = 2
tss_c4 = root1 * "Segmentation_type6/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, channel; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)



ROOT1 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Ifnb1/Exp1/"
root1 = ROOT1
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Ifnb1_20210125.csv", DataFrame)
expn = 1
gene = "Ifnb1forL2"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Ifnb1L2")
gene = "L2"
typeprobe = "type6"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Ifnb1L2")




function get_image_patterns2(imagefolder)
    files = readdir(imagefolder)
    bool = [length(split(f, r"_C[\d][\.].|_C[\d]_")) > 1 for f in files]
    files = files[bool]
    return unique([split(f, r"_C[\d][\.].|_C[\d]_")[1] for f in files])
end

function TSS_raw_quant2(t2, tss_folder, image_folder, n; xy = 0.189, zx = 0.5)
    
    images_pat = get_image_patterns2(t2)
    
    p = Progress(length(images_pat), 1)

    dfs = []

    for a in 1:length(images_pat)
        next!(p)
        

        pat = images_pat[a]

        # Get all the TSS in the images
        a = TSSs.find_outline(tss_folder, pat)
        tss = TSSs.cell_tss_dict(a)
        tss = TSSs.meassure_tss(tss, imagesfolder, pat, n, xy = xy, zx = zx)
    
        df = TSSs.analysis_singlechannel(tss; image = pat,  tss_name = :TSS2)
        push!(dfs, df)

    end
    
    d = TSSs.join_in_all_common_columns(dfs...)
    return d
    
end

ROOT2 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Ifnb1/Exp2/"

root1 = ROOT2
expn = 2
gene = "Ifnb1forL2"
typeprobe = "type4"
channel= 3
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = get_image_patterns2(tss_c4)
d = TSS_raw_quant2(tss_c4, tss_c4, imagesfolder, channel; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)


gene = "L2"
typeprobe = "type6"
channel = 2
tss_c4 = root1 * "Segmentation_type6/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = get_image_patterns2(tss_c4)
d = TSS_raw_quant2(tss_c4, tss_c4, imagesfolder, channel; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

function FQfiles.fix_image_name(df)
    # Get rid of the channel specific identifiers
    df[!,:Image] = [split(i, r"_C[\d][\.].|_C[\d]_")[1] for i in df[!,:Image]]
    df
end

ROOT2 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Ifnb1/Exp2/"
root1 = ROOT2
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Ifnb1_20210205.csv", DataFrame)
expn = 2
gene = "Ifnb1forL2"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Ifnb1L2")
gene = "L2"
typeprobe = "type6"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Ifnb1L2")