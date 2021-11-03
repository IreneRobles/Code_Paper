if !in("Il12b_original", readdir())
    mkdir("Il12b_original")
end

function normalise_avgdot(df, genename, expn, typeprobe, ibidislide; root = pwd(), genefolder = "Gene")
    dot1_r1 = TSSs.int_brightest_pixel(TSSs.read_tiff_as_gray(root * typeprobe *"/_mRNA_AVG_ns.tif").*255; radious = 1)
    dot1_r2 = TSSs.int_brightest_pixel(TSSs.read_tiff_as_gray(root * typeprobe * "/_mRNA_AVG_ns.tif").*255; radious = 2)
    dot1_r3 = TSSs.int_brightest_pixel(TSSs.read_tiff_as_gray(root * typeprobe * "/_mRNA_AVG_ns.tif").*255; radious = 3)

    df[!,:TSS1_r1] = df[!,:locus1_int1_TSS2] ./ dot1_r1
    df[!,:TSS1_r2] = df[!,:locus1_int2_TSS2] ./ dot1_r2
    df[!,:TSS1_r3] = df[!,:locus1_int3_TSS2] ./ dot1_r3
    df[!,:TSS2_r1] = df[!,:locus2_int1_TSS2] ./ dot1_r1
    df[!,:TSS2_r2] = df[!,:locus2_int2_TSS2] ./ dot1_r2
    df[!,:TSS2_r3] = df[!,:locus2_int3_TSS2] ./ dot1_r3
    
    CSV.write("TSS_avgdot/"*genename*"_exp"*"$expn"*".csv", df)
    ibidislide[!,:Well] = [split(ii, " (")[1] for ii in ibidislide[!,:Well]]
    exp1 = CSV.read("TSS_avgdot/"*gene*"_exp"*"$expn"*".csv", DataFrame)
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_TSS.csv", exp1)
    fq_dir1 = root1*typeprobe*"/"
    fq1 = FQ_summary_MATURE(fq_dir1)
    CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ.csv", fq1)
end

function timepoint(ii)
    if occursin("U 0", ii) 
        return "0"
    elseif occursin("U 1", ii) 
        return "60"
    elseif occursin("U 2", ii) 
        return "90"
    elseif occursin("U 3", ii) || occursin("FL120", ii)
        return "120"
    else 
        "120"
    end
end


function genotype(ii)
    if occursin("V 0", ii) 
        return "WT"
    elseif occursin("V 1", ii) 
        return "Rad21KO" 
     elseif occursin("FL120", ii)
        return "Rad21KO"
    else 
        "unknown"
    end
end

ROOT1 = "/Volumes/lymphdev\$/IreneR/Confocal/il12b_intron_exon/Exp1/"
ROOT2 = "/Volumes/lymphdev\$/IreneR/Confocal/il12b_intron_exon/Exp2/"
ROOT3 = "/Volumes/lymphdev\$/IreneR/Confocal/il12b_intron_exon/Exp4/"
ROOT4 = "/Volumes/lymphdev\$/IreneR/Confocal/il12b_intron_exon/Exp6/"


root1 = ROOT1
expn = 1
gene = "Il12b_intron_original"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)


root1 = ROOT2
expn = 2
gene = "Il12b_intron_original"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT3
expn = 3
gene = "Il12b_intron_original"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

root1 = ROOT4
expn = 4
gene = "Il12b_intron_original"
tss_c4 = root1 * "Segmentation_type4/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 3; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp"*"$expn"*".csv", d)

folder = "Il12b_original"

root1 = ROOT1
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Il12boriginal_1.csv", DataFrame)
expn = 1
gene = "Il12b_intron_original"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = folder)

root1 = ROOT2
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Il12boriginal_2.csv", DataFrame)
expn = 2
gene = "Il12b_intron_original"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = folder)

root1 = ROOT3
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Il12boriginal_3.csv", DataFrame)
expn = 3
gene = "Il12b_intron_original"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = folder)

root1 = ROOT4
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Il12boriginal_4.csv", DataFrame)
expn = 4
gene = "Il12b_intron_original"
typeprobe = "type4"
tss1 = CSV.read("TSS_raw/"*gene*"_exp"*"$expn"*".csv", DataFrame)
normalise_avgdot(tss1, gene, expn, typeprobe, exp; root = root1, genefolder = folder)

genefolder = "Il12b_original"
root1 = ROOT1
expn = 1
gene = "Il12b_intron_original"
typeprobe = "type6"
fq_dir1 = root1*typeprobe*"/"
fq1 = FQ_summary_MATURE(fq_dir1)
CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ_mature.csv", fq1)

root1 = ROOT2
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Il12boriginal_2.csv", DataFrame)
expn = 2
gene = "Il12b_intron_original"
typeprobe = "type6"
fq_dir1 = root1*typeprobe*"/"
fq1 = FQ_summary_MATURE(fq_dir1)
CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ_mature.csv", fq1)

root1 = ROOT3
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Il12boriginal_3.csv", DataFrame)
expn = 3
gene = "Il12b_intron_original"
typeprobe = "type6"
fq_dir1 = root1*typeprobe*"/"
fq1 = FQ_summary_MATURE(fq_dir1)
CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ_mature.csv", fq1)

root1 = ROOT4
exp = CSV.read("IbidiChambers/IbidiChamberSlide_Il12boriginal_4.csv", DataFrame)
expn = 4
gene = "Il12b_intron_original"
typeprobe = "type6"
fq_dir1 = root1*typeprobe*"/"
fq1 = FQ_summary_MATURE(fq_dir1)
CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ_mature.csv", fq1)

folder = "Il12b_original"

expn = 1
gene = "Il12b_intron_original"
tss1 = CSV.read(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ_mature.csv", DataFrames.DataFrame)
tss1[!,:Timepoint] = [ timepoint(ii) for ii in tss1[!,:Image]]
tss1[!,:Genotype] = [ genotype(ii) for ii in tss1[!,:Image]]
CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_CP.csv", tss1[!,["Genotype","Timepoint","Image_Cell"]])

expn = 2
gene = "Il12b_intron_original"
tss1 = CSV.read(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ_mature.csv", DataFrames.DataFrame)
tss1[!,:Timepoint] = [ timepoint(ii) for ii in tss1[!,:Image]]
tss1[!,:Genotype] = [ genotype(ii) for ii in tss1[!,:Image]]
CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_CP.csv", tss1[!,["Genotype","Timepoint","Image_Cell"]])

expn = 3
gene = "Il12b_intron_original"
tss1 = CSV.read(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ_mature.csv", DataFrames.DataFrame)
tss1[!,:Timepoint] = [ timepoint(ii) for ii in tss1[!,:Image]]
tss1[!,:Genotype] = [ genotype(ii) for ii in tss1[!,:Image]]
CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_CP.csv", tss1[!,["Genotype","Timepoint","Image_Cell"]])

expn = 4
gene = "Il12b_intron_original"
tss1 = CSV.read(genefolder*"/"*gene*"_exp"*"$expn"*"_FQ_mature.csv", DataFrames.DataFrame)
tss1[!,:Timepoint] = [ timepoint(ii) for ii in tss1[!,:Image]]
tss1[!,:Genotype] = [ genotype(ii) for ii in tss1[!,:Image]]
CSV.write(genefolder*"/"*gene*"_exp"*"$expn"*"_CP.csv", tss1[!,["Genotype","Timepoint","Image_Cell"]])


