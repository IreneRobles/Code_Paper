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

folder = "Prdm1"
ROOT1 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Prdm1/Exp1/"
ROOT2 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Prdm1/Exp2/"

root1 = ROOT1

expn = 1
gene = "Prdm1_enh"
typeprobe = "type6"
tss_c4 = root1 * "Segmentation_type6/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT2
expn = 2
gene = "Prdm1_enh"
typeprobe = "type6"
tss_c4 = root1 * "Segmentation_type6/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT1

expn = 1
gene = "Prdm1_intron"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT2
expn = 2
gene = "Prdm1_intron"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT1
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Prdm1_20210125.csv", DataFrame)
expn = 1
gene = "Prdm1_enh"
typeprobe = "type6"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = folder)

root1 = ROOT2
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Prdm1_20210205.csv", DataFrame)
expn = 2
gene = "Prdm1_enh"
typeprobe = "type6"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = folder)

root1 = ROOT1
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Prdm1_20210125.csv", DataFrame)
expn = 1
gene = "Prdm1_intron"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = folder)

root1 = ROOT2
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Prdm1_20210205.csv", DataFrame)
expn = 2
gene = "Prdm1_intron"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = folder)