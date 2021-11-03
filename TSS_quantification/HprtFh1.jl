  if !in("Fh1", readdir())
    mkdir("Fh1")
    end
  if !in("Hprt", readdir())
    mkdir("Hprt")
    end


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

function timepoint(ii)
    if occursin("U 0", ii) 
        return "0"
    elseif occursin("U 1", ii) 
        return "0i"
    elseif occursin("U 2", ii) 
        return "8"
    elseif occursin("U 3", ii) 
        return "8i"
    else 
        "0"
    end
end


function genotype(ii)
    if occursin("V 0", ii) 
        return "WT"
    elseif occursin("V 1", ii) 
        return "Rad21KO"
    end
end

# Define where are the outlines

root = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp1"
tss_c2 = root * "/Segmentation_type6/_FQ_outline/_TS_detect"
tss_c3 = root * "/Segmentation_type4/_FQ_outline/_TS_detect"
tss_c4 = root * "/Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c2)

if !in("TSS_raw", readdir())
    mkdir("TSS_raw")
end

d = TSSs.TSS_raw_quant(tss_c2, tss_c2, imagesfolder, 2)

d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]
CSV.write("TSS_raw/Hprt_info.csv", d)

d = TSSs.TSS_raw_quant(tss_c2, tss_c3, imagesfolder, 3)
d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]
CSV.write("TSS_raw/Fh1_info.csv", d)

d = TSSs.TSS_raw_quant(tss_c2, tss_c4, imagesfolder, 4)
d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]
CSV.write("TSS_raw/Actb_info.csv", d)


root = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp2"
tss_c2 = root * "/Segmentation_type6/_FQ_outline/_TS_detect"
tss_c3 = root * "/Segmentation_type4/_FQ_outline/_TS_detect"
tss_c4 = root * "/Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c3)

d = TSSs.TSS_raw_quant(tss_c2, tss_c2, imagesfolder, 2)

d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]
CSV.write("TSS_raw/Hprt_info2.csv", d)


d = TSSs.TSS_raw_quant(tss_c2, tss_c3, imagesfolder, 3)
d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]
CSV.write("TSS_raw/Fh1_info2.csv", d)


d = TSSs.TSS_raw_quant(tss_c2, tss_c4, imagesfolder, 4)
d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]
CSV.write("TSS_raw/Actb_info2.csv", d)


root = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp3"
tss_c2 = root * "/Segmentation_type6/_FQ_outline/_TS_detect"
tss_c3 = root * "/Segmentation_type4/_FQ_outline/_TS_detect"
tss_c4 = root * "/Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c3)

if !in("TSS_raw", readdir())
    mkdir("TSS_raw")
end

d = TSSs.TSS_raw_quant(tss_c2, tss_c2, imagesfolder, 2)
d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]
CSV.write("TSS_raw/Hprt_info3.csv", d)

d = TSSs.TSS_raw_quant(tss_c2, tss_c3, imagesfolder, 3)
d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]

CSV.write("TSS_raw/Fh1_info3.csv", d)

d = TSSs.TSS_raw_quant(tss_c2, tss_c4, imagesfolder, 4)
d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]

CSV.write("TSS_raw/Actb_info3.csv", d)


# Define where are the outlines

# Define where are the outlines

root = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp6"
tss_c2 = root * "/Segmentation_type6/_FQ_outline/_TS_detect"
tss_c3 = root * "/Segmentation_type4/_FQ_outline/_TS_detect"
tss_c4 = root * "/Segmentation_type1/_FQ_outline/_TS_detect"
imagesfolder = root * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c3)

d = TSSs.TSS_raw_quant(tss_c2, tss_c2, imagesfolder, 2)

d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]
CSV.write("TSS_raw/Hprt_info6.csv", d)
d = TSSs.TSS_raw_quant(tss_c2, tss_c3, imagesfolder, 3)


d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]


CSV.write("TSS_raw/Fh1_info6.csv", d)

function TSSs.TSS_raw_quant(t2, tss_folder, image_folder, n; xy = 0.189, zx = 0.5)
    
    images_pat = intersect(TSSs.get_image_patterns(t2), TSSs.get_image_patterns(tss_folder))
    
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

d = TSSs.TSS_raw_quant(tss_c2, tss_c4, imagesfolder, 4)


d[!,:Timepoint] = [timepoint(ii) for ii in d[!,:Image]]
d[!,:Genotype] = [genotype(ii) for ii in d[!,:Image]]

CSV.write("TSS_raw/Actb_info6.csv", d)


orgdatafold = ""

f_hprt = orgdatafold*"Hprt"
f_fh1 = orgdatafold*"Fh1"


root1 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp1"
root2 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp2"
root3 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp3"
root4 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp6"
if !in("TSS_avgdot", readdir())
    mkdir("TSS_avgdot")
end

if !in("TSS_avgdot", readdir())
    mkdir("TSS_avgdot")
end

gene = "Hprt"

tss1 = CSV.read("TSS_raw/"*gene*"_info.csv", DataFrame)
tss2 = CSV.read("TSS_raw/"*gene*"_info2.csv", DataFrame)
tss3 = CSV.read("TSS_raw/"*gene*"_info3.csv", DataFrame)
tss4 = CSV.read("TSS_raw/"*gene*"_info6.csv", DataFrame)


dot1_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 1)
dot2_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 1)
dot3_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 1)
dot4_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root4 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 1)

dot1_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 2)
dot2_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 2)
dot3_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 2)
dot4_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root4 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 2)

dot1_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 3)
dot2_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 3)
dot3_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 3)
dot4_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root4 * "/type6/_mRNA_AVG_ns.tif").*nor; radious = 3)

tss1[!,:TSS1_r1] = tss1[!,:locus1_int1_TSS2] ./ dot2_r1
tss1[!,:TSS1_r2] = tss1[!,:locus1_int2_TSS2] ./ dot2_r2
tss1[!,:TSS1_r3] = tss1[!,:locus1_int3_TSS2] ./ dot2_r3


tss2[!,:TSS1_r1] = tss2[!,:locus1_int1_TSS2] ./ dot2_r1
tss2[!,:TSS1_r2] = tss2[!,:locus1_int2_TSS2] ./ dot2_r2
tss2[!,:TSS1_r3] = tss2[!,:locus1_int2_TSS2] ./ dot2_r3

tss3[!,:TSS1_r1] = tss3[!,:locus1_int1_TSS2] ./ dot3_r1
tss3[!,:TSS1_r2] = tss3[!,:locus1_int2_TSS2] ./ dot3_r2
tss3[!,:TSS1_r3] = tss3[!,:locus1_int3_TSS2] ./ dot3_r3


tss4[!,:TSS1_r1] = tss4[!,:locus1_int1_TSS2] ./ dot4_r1
tss4[!,:TSS1_r2] = tss4[!,:locus1_int2_TSS2] ./ dot4_r2
tss4[!,:TSS1_r3] = tss4[!,:locus1_int3_TSS2] ./ dot4_r3


tss1[!,:TSS2_r1] = tss1[!,:locus2_int1_TSS2] ./ dot2_r1
tss1[!,:TSS2_r2] = tss1[!,:locus2_int2_TSS2] ./ dot2_r2
tss1[!,:TSS2_r3] = tss1[!,:locus2_int3_TSS2] ./ dot2_r3

tss2[!,:TSS2_r1] = tss2[!,:locus2_int1_TSS2] ./ dot2_r1
tss2[!,:TSS2_r2] = tss2[!,:locus2_int2_TSS2] ./ dot2_r2
tss2[!,:TSS2_r3] = tss2[!,:locus2_int3_TSS2] ./ dot2_r3

tss3[!,:TSS2_r1] = tss3[!,:locus2_int1_TSS2] ./ dot3_r1
tss3[!,:TSS2_r2] = tss3[!,:locus2_int2_TSS2] ./ dot3_r2
tss3[!,:TSS2_r3] = tss3[!,:locus2_int3_TSS2] ./ dot3_r3


tss4[!,:TSS2_r1] = tss4[!,:locus2_int1_TSS2] ./ dot4_r1
tss4[!,:TSS2_r2] = tss4[!,:locus2_int2_TSS2] ./ dot4_r2
tss4[!,:TSS2_r3] = tss4[!,:locus2_int3_TSS2] ./ dot4_r3

CSV.write(orgdatafold*"TSS_avgdot/"*gene*"_exp1.csv", tss1)
CSV.write(orgdatafold*"TSS_avgdot/"*gene*"_exp2.csv", tss2)
CSV.write(orgdatafold*"TSS_avgdot/"*gene*"_exp3.csv", tss3)
CSV.write(orgdatafold*"TSS_avgdot/"*gene*"_exp4.csv", tss4)


root1 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp1"
root2 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp2"
root3 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp3"
root4 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp6"
if !in("TSS_avgdot", readdir())
    mkdir("TSS_avgdot")
end

if !in("TSS_avgdot", readdir())
    mkdir("TSS_avgdot")
end

gene = "Fh1"

tss1 = CSV.read("TSS_raw/"*gene*"_info.csv", DataFrame)
tss2 = CSV.read("TSS_raw/"*gene*"_info2.csv", DataFrame)
tss3 = CSV.read("TSS_raw/"*gene*"_info3.csv", DataFrame)
tss4 = CSV.read("TSS_raw/"*gene*"_info6.csv", DataFrame)


dot1_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 1)
dot2_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 1)
dot3_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 1)
dot4_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root4 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 1)

dot1_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 2)
dot2_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 2)
dot3_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 2)
dot4_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root4 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 2)

dot1_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 3)
dot2_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 3)
dot3_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 3)
dot4_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root4 * "/type4/_mRNA_AVG_ns.tif").*nor; radious = 3)

tss1[!,:TSS1_r1] = tss1[!,:locus1_int1_TSS2] ./ dot2_r1
tss1[!,:TSS1_r2] = tss1[!,:locus1_int2_TSS2] ./ dot2_r2
tss1[!,:TSS1_r3] = tss1[!,:locus1_int3_TSS2] ./ dot2_r3


tss2[!,:TSS1_r1] = tss2[!,:locus1_int1_TSS2] ./ dot2_r1
tss2[!,:TSS1_r2] = tss2[!,:locus1_int2_TSS2] ./ dot2_r2
tss2[!,:TSS1_r3] = tss2[!,:locus1_int2_TSS2] ./ dot2_r3

tss3[!,:TSS1_r1] = tss3[!,:locus1_int1_TSS2] ./ dot3_r1
tss3[!,:TSS1_r2] = tss3[!,:locus1_int2_TSS2] ./ dot3_r2
tss3[!,:TSS1_r3] = tss3[!,:locus1_int3_TSS2] ./ dot3_r3


tss4[!,:TSS1_r1] = tss4[!,:locus1_int1_TSS2] ./ dot4_r1
tss4[!,:TSS1_r2] = tss4[!,:locus1_int2_TSS2] ./ dot4_r2
tss4[!,:TSS1_r3] = tss4[!,:locus1_int3_TSS2] ./ dot4_r3


tss1[!,:TSS2_r1] = tss1[!,:locus2_int1_TSS2] ./ dot2_r1
tss1[!,:TSS2_r2] = tss1[!,:locus2_int2_TSS2] ./ dot2_r2
tss1[!,:TSS2_r3] = tss1[!,:locus2_int3_TSS2] ./ dot2_r3

tss2[!,:TSS2_r1] = tss2[!,:locus2_int1_TSS2] ./ dot2_r1
tss2[!,:TSS2_r2] = tss2[!,:locus2_int2_TSS2] ./ dot2_r2
tss2[!,:TSS2_r3] = tss2[!,:locus2_int3_TSS2] ./ dot2_r3

tss3[!,:TSS2_r1] = tss3[!,:locus2_int1_TSS2] ./ dot3_r1
tss3[!,:TSS2_r2] = tss3[!,:locus2_int2_TSS2] ./ dot3_r2
tss3[!,:TSS2_r3] = tss3[!,:locus2_int3_TSS2] ./ dot3_r3


tss4[!,:TSS2_r1] = tss4[!,:locus2_int1_TSS2] ./ dot4_r1
tss4[!,:TSS2_r2] = tss4[!,:locus2_int2_TSS2] ./ dot4_r2
tss4[!,:TSS2_r3] = tss4[!,:locus2_int3_TSS2] ./ dot4_r3

CSV.write(orgdatafold*"TSS_avgdot/"*gene*"_exp1.csv", tss1)
CSV.write(orgdatafold*"TSS_avgdot/"*gene*"_exp2.csv", tss2)
CSV.write(orgdatafold*"TSS_avgdot/"*gene*"_exp3.csv", tss3)
CSV.write(orgdatafold*"TSS_avgdot/"*gene*"_exp4.csv", tss4)



function get_cells(root)
    tss_c2 = root * "/Segmentation_type6/_FQ_outline/_TS_detect"
    tss_c3 = root * "/Segmentation_type4/_FQ_outline/_TS_detect"
    tss_c4 = root * "/Segmentation_type1/_FQ_outline/_TS_detect"

    imagesfolder = root * "/tiff3D"

    cp_folder = root *"/CP_results"
    images = CSV.read(cp_folder * "/MyExpt_Image.csv", DataFrame)
    images = CSV.read(cp_folder * "/MyExpt_Image.csv", DataFrame)
    cells = CSV.read(cp_folder * "/MyExpt_TrueCells.csv", DataFrame);
    

    images[!,:Image] = [split(ii, "_C")[1] for ii in images[!,:FileName_DAPI]]
    images[!,:Genotype] = [genotype(ii) for ii in images[!,:FileName_DAPI]]
    images[!,:Timepoint] = [timepoint(ii) for ii in images[!,:FileName_DAPI]]
    images[!,:Sample] = images[!,:Genotype] .* "_" .* string.(images[!,:Timepoint]) 
    im = images[!, [:Sample, :ImageNumber, :Genotype, :Timepoint, :Image]]

    c1 = innerjoin(cells, im, on = :ImageNumber);

    c1[!,:Cell_Image] = "Cell_CP_" .* string.(c1[!,:ObjectNumber]) .* "__" .* c1[!,:Image]
    c1[!,:Cell] = "Cell_CP_" .* string.(c1[!,:ObjectNumber])
    
        c1 = column_fusion(c1, :Image, :Cell)

    return c1
end

function get_TSS_data(gene, cells, rep)
    suff = if rep == 1 "" else string(rep) end
    file = orgdatafold*"TSS_avgdot/"*gene*"_exp"*suff*".csv"
    dat5 = CSV.read(file, DataFrame)
    dat5[!,:Cell_Image] = dat5[!,:Cell].* "__" .* dat5[!,:Image]
    dat5 = innerjoin(cells[!, [:Cell_Image]], dat5, on = :Cell_Image)
    dat5[!,:Genotype] = [genotype(ii) for ii in dat5[!,:Image]]
    dat5[!,:Timepoint] = [timepoint(ii) for ii in dat5[!,:Image]]
    dat5[!,:Sample] = dat5[!,:Genotype] .* "_" .* string.(dat5[!,:Timepoint])
    dat5 = column_fusion(dat5, :Image, :Cell)

    return dat5
    
end

gene = "Hprt"

exp1 = get_TSS_data(gene, get_cells(root1), "1")
exp2 = get_TSS_data(gene, get_cells(root2), "2")
# WT in e2 is weird 
exp2 = exp2[exp2[!,:Sample].!= "WT_0", :]

exp3 = get_TSS_data(gene, get_cells(root3), "3")
exp4 = get_TSS_data(gene, get_cells(root4), "4")
exp4 = exp4[exp4[!,:Sample].!= "Rad21KO_8i", :]


CSV.write(orgdatafold*gene*"/"*gene*"_exp1_TSS.csv", exp1)
CSV.write(orgdatafold*gene*"/"*gene*"_exp1_CP.csv", get_cells(root1))

CSV.write(orgdatafold*gene*"/"*gene*"_exp2_TSS.csv", exp2)
CSV.write(orgdatafold*gene*"/"*gene*"_exp2_CP.csv", get_cells(root2))

CSV.write(orgdatafold*gene*"/"*gene*"_exp3_TSS.csv", exp3)
CSV.write(orgdatafold*gene*"/"*gene*"_exp3_CP.csv", get_cells(root3))

CSV.write(orgdatafold*gene*"/"*gene*"_exp4_TSS.csv", exp4)
CSV.write(orgdatafold*gene*"/"*gene*"_exp4_CP.csv", get_cells(root4))


gene = "Fh1"

exp1 = get_TSS_data(gene, get_cells(root1), "1")
exp2 = get_TSS_data(gene, get_cells(root2), "2")
# WT in e2 is weird 
exp2 = exp2[exp2[!,:Sample].!= "WT_0", :]


exp3 = get_TSS_data(gene, get_cells(root3), "3")
exp4 = get_TSS_data(gene, get_cells(root4), "4")
exp4 = exp4[exp4[!,:Sample].!= "Rad21KO_8i", :]


CSV.write(orgdatafold*gene*"/"*gene*"_exp1_TSS.csv", exp1)
CSV.write(orgdatafold*gene*"/"*gene*"_exp1_CP.csv", get_cells(root1))

CSV.write(orgdatafold*gene*"/"*gene*"_exp2_TSS.csv", exp2)
CSV.write(orgdatafold*gene*"/"*gene*"_exp2_CP.csv", get_cells(root2))

CSV.write(orgdatafold*gene*"/"*gene*"_exp3_TSS.csv", exp3)
CSV.write(orgdatafold*gene*"/"*gene*"_exp3_CP.csv", get_cells(root3))

CSV.write(orgdatafold*gene*"/"*gene*"_exp4_TSS.csv", exp4)
CSV.write(orgdatafold*gene*"/"*gene*"_exp4_CP.csv", get_cells(root4))


root1 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp1"
root2 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp2"
root3 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp3"
root4 = "/Volumes/lymphdev\$/IreneR/Confocal/HprtFh1Actb/Exp6"


gene = "Fh1"
fq1 = FQ_summary_MATURE(root1*"/type4/")
CSV.write(orgdatafold*gene*"/"*gene*"_exp1_FQ.csv", fq1)
fq1 = FQ_summary_MATURE(root2*"/type4/")
CSV.write(orgdatafold*gene*"/"*gene*"_exp2_FQ.csv", fq1)
fq1 = FQ_summary_MATURE(root3*"/type4/")
CSV.write(orgdatafold*gene*"/"*gene*"_exp3_FQ.csv", fq1)
fq1 = FQ_summary_MATURE(root4*"/type4/")
CSV.write(orgdatafold*gene*"/"*gene*"_exp4_FQ.csv", fq1)

gene = "Hprt"
fq1 = FQ_summary_MATURE(root1*"/type6/")
CSV.write(orgdatafold*gene*"/"*gene*"_exp1_FQ.csv", fq1)
fq1 = FQ_summary_MATURE(root2*"/type6/")
CSV.write(orgdatafold*gene*"/"*gene*"_exp2_FQ.csv", fq1)
fq1 = FQ_summary_MATURE(root3*"/type6/")
CSV.write(orgdatafold*gene*"/"*gene*"_exp3_FQ.csv", fq1)
fq1 = FQ_summary_MATURE(root4*"/type6/")
CSV.write(orgdatafold*gene*"/"*gene*"_exp4_FQ.csv", fq1)