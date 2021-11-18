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

ROOT1 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Il12bHSS1timecourse/Rep1_2/"
ROOT2 = "/Volumes/lymphdev\$/IreneR/Confocal/EnhancerInducibleGenes/Il12bHSS1timecourse/Rep3_4/"

root1 = ROOT1
expn = "1"

gene = "HSS1_timecourse"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

gene = "Il12b_intron_timecourse"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

root1 = ROOT2
expn = "2"

gene = "HSS1_timecourse"
typeprobe = "type1"
tss_c4 = root1 * "Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

gene = "Il12b_intron_timecourse"
typeprobe = "type4"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)

folder = "Il12bHSS1BD"
root1 = ROOT1
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_timecourse_exp1.csv", DataFrame)
expn = 1
gene = "HSS1_timecourse"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")

gene = "Il12b_intron_timecourse"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")

root1 = ROOT2
exp = CSV.read("IbidiChambers/IbidiChamber_Il12bHSS1_timecourse_exp2.csv", DataFrame)
expn = 2
gene = "HSS1_timecourse"
typeprobe = "type1"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")

gene = "Il12b_intron_timecourse"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = "Il12bHSS1")
