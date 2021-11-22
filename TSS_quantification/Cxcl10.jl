if !in("Cxcl10", readdir())
    mkdir("Cxcl10")
end


function CellInfo(dir, exp_df)
    image_file = get_files_ending_with(dir, "Image.csv")
     cell_file = get_files_ending_with(dir, "TrueCells.csv")

    if .&(length(cell_file) == 1,  length(image_file) == 1) # Make sure that there is only one file with the image information
        image_file = normpath(dir, image_file[1])
        cell_file = normpath(dir, cell_file[1])
    else
        error("Image.csv or Cell file not found")
    end
    
    ima = DataFrames.DataFrame(CSV.read(image_file))[:, [:FileName_DAPI, :ImageNumber]]
    ima[:Image] = [split(split(ii, "_C1")[1], "_MAX")[1] for ii in ima[:FileName_DAPI]]
    
    ima[:WELL] = [split(ii, " ") for ii in ima[:Image]]
    ima[:WELL] = [1 for ii in ima[:WELL]]
    
     ima[:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in ima[:FileName_DAPI]]
    
    exp_df[:Well] = [split(ii, " (")[1] for ii in exp_df[:Well]]

    im_ = join(exp_df, ima, on = :Well, makeunique=true)
        

    
   
     cells = DataFrames.DataFrame(CSV.read(cell_file))
    cells = join(im_, cells, on = :ImageNumber)
    
    cells[:Image_Cell] = [cells[ii, :Image]*"__Cell_CP_"*string(cells[ii, :ObjectNumber]) for ii in 1:nrow(cells)]
    return dropmissing(cells)
    
end


function NuInfo(dir, exp_df)
    image_file = get_files_ending_with(dir, "Image.csv")
     cell_file = get_files_ending_with(dir, "TrueNuclei.csv")

    if .&(length(cell_file) == 1,  length(image_file) == 1) # Make sure that there is only one file with the image information
        image_file = normpath(dir, image_file[1])
        cell_file = normpath(dir, cell_file[1])
    else
        error("Image.csv or Cell file not found")
    end
    
    ima = DataFrames.DataFrame(CSV.read(image_file))[:, [:FileName_DAPI, :ImageNumber]]
    ima[:Image] = [split(split(ii, "_C1")[1], "_MAX")[1] for ii in ima[:FileName_DAPI]]
    
    ima[:WELL] = [split(ii, " ") for ii in ima[:Image]]
    ima[:WELL] = [1 for ii in ima[:WELL]]
    
     ima[:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in ima[:FileName_DAPI]]
    
    exp_df[:Well] = [split(ii, " (")[1] for ii in exp_df[:Well]]

    im_ = join(exp_df, ima, on = :Well, makeunique=true)
        

    
   
     cells = DataFrames.DataFrame(CSV.read(cell_file))
    cells = join(im_, cells, on = :ImageNumber)
    
    cells[:Image_Cell] = [cells[ii, :Image]*"__Cell_CP_"*string(cells[ii, :ObjectNumber]) for ii in 1:nrow(cells)]
    return dropmissing(cells)
    
end



function CellInfo1(dir)
    image_file = get_files_ending_with(dir, "Image.csv")
    cell_file = get_files_ending_with(dir, "TrueCells.csv")

    if .&(length(cell_file) == 1,  length(image_file) == 1) # Make sure that there is only one file with the image information
        image_file = string(dir, image_file[1])
        cell_file = string(dir, cell_file[1])
    else
        error("Image.csv or Cell file not found")
    end
    
     ima = CSV.read(image_file, DataFrames.DataFrame)[!, [:FileName_DAPI, :ImageNumber]]
    ima[!, :Image] = [split(split(ii, "_C1.")[1], "MAX_")[2] for ii in ima[!,:FileName_DAPI]]
    info = [split(ii, "_") for ii in Array(ima[!, :Image])]    
    ima[!, :Timepoint] = extract_metainformation(info, 
        ["V 0", "V 1"], 
        [0, 120])  
    ima[!, :Plate] = [ii[3] for ii in info]

    # Describe samples present in control plate
    ima_controls = ima[ima[!,:Plate].== "Controls", :]
    info_controls  = [split(ii, "_") for ii in Array(ima_controls[!,:Image])]   
    ima_controls[!, :Genotype] = extract_metainformation(info_controls, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Ifnb1KO", "WT","Rad21KO", "Ifnb1KO"])
    ima_controls[!, :Probe] = extract_metainformation(info_controls, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Cxcl10", "Cxcl10","Cxcl10", "Mx1"])
    
     # Describe samples present in mix plate
    ima_mix = ima[ima[!, :Plate].== "Mix", :]
    info_mix  = [split(ii, "_") for ii in Array(ima_mix[!, :Image])]   
    ima_mix[!, :Genotype] = extract_metainformation(info_mix, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Ifnb1KO+WT", "Ifnb1KO+Rad21KO","Ifnb1KO+WT", "Ifnb1KO+Rad21KO"])
    ima_mix[!, :Probe] = extract_metainformation(info_mix, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Cxcl10", "Cxcl10","Mx1", "Mx1"])
        
    imdata = join_in_all_common_columns(ima_mix, ima_controls)
    
     imdata[!, :Sample] = [imdata[!,:Genotype][ii]*"_"*string(imdata[!, :Timepoint][ii])*"_"*string(imdata[!, :Probe][ii]) for ii in 1:nrow(imdata)]

    
    csv = CSV.read(cell_file, DataFrames.DataFrame)
    csv[!, :Rep] = [1 for i in 1:nrow(csv)]
    csv = innerjoin(csv, imdata, on = :ImageNumber)
    
    csv[!, :Image_Cell] = [csv[ii, :Image]*"__Cell_CP_"*string(csv[ii, :ObjectNumber]) for ii in 1:nrow(csv)]
    return csv
end

function CellInfo(dir)
    image_file = get_files_ending_with(dir, "Image.csv")
    cell_file = get_files_ending_with(dir, "TrueCells.csv")

    if .&(length(cell_file) == 1,  length(image_file) == 1) # Make sure that there is only one file with the image information
        image_file = string(dir, image_file[1])
        cell_file = string(dir, cell_file[1])
    else
        error("Image.csv or Cell file not found")
    end
    
     ima = CSV.read(image_file, DataFrame)[!, [:FileName_DAPI, :ImageNumber]]
    ima[!, :Image] = [split(split(ii, "_C1.")[1], "MAX_")[2] for ii in ima[!,:FileName_DAPI]]
    info = [split(ii, "_") for ii in Array(ima[!, :Image])]    
    ima[!, :Timepoint] = extract_metainformation(info, 
        ["V 0", "V 1"], 
        [0, 120])  
    ima[!, :Plate] = [ii[3] for ii in info]

    # Describe samples present in control plate
    ima_controls = ima[ima[!,:Plate].== "Controls", :]
    info_controls  = [split(ii, "_") for ii in Array(ima_controls[!,:Image])]   
    ima_controls[!, :Genotype] = extract_metainformation(info_controls, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Ifnb1KO", "Rad21KO","WT", "Ifnb1KO"])
    ima_controls[!, :Probe] = extract_metainformation(info_controls, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Cxcl10", "Cxcl10","Cxcl10", "Mx1"])
    
     # Describe samples present in mix plate
    ima_mix = ima[ima[!, :Plate].== "Mix", :]
    info_mix  = [split(ii, "_") for ii in Array(ima_mix[!, :Image])]   
    ima_mix[!, :Genotype] = extract_metainformation(info_mix, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Ifnb1KO+WT", "Ifnb1KO+Rad21KO","Ifnb1KO+WT", "Ifnb1KO+Rad21KO"])
    ima_mix[!, :Probe] = extract_metainformation(info_mix, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Cxcl10", "Cxcl10","Mx1", "Mx1"])
        
    imdata = join_in_all_common_columns(ima_mix, ima_controls)
    
     imdata[!, :Sample] = [imdata[!,:Genotype][ii]*"_"*string(imdata[!, :Timepoint][ii])*"_"*string(imdata[!, :Probe][ii]) for ii in 1:nrow(imdata)]

    
    csv = CSV.read(cell_file, DataFrames.DataFrame)
    csv[!, :Rep] = [1 for i in 1:nrow(csv)]
    csv = innerjoin(csv, imdata, on = :ImageNumber)
    
    csv[!, :Image_Cell] = [csv[ii, :Image]*"__Cell_CP_"*string(csv[ii, :ObjectNumber]) for ii in 1:nrow(csv)]
    return csv
end


function CellInfo2(dir)
    image_file = get_files_ending_with(dir, "Image.csv")
     cell_file = get_files_ending_with(dir, "TrueCells.csv")

    if .&(length(cell_file) == 1,  length(image_file) == 1) # Make sure that there is only one file with the image information
        image_file = string(dir, image_file[1])
        cell_file = string(dir, cell_file[1])
    else
        error("Image.csv or Cell file not found")
    end
    
     ima = CSV.read(image_file, DataFrames.DataFrame)[:, [:FileName_DAPI, :ImageNumber]]
    ima[!, :Image] = [split(split(ii, "_C1.")[1], "MAX_")[2] for ii in ima[!,:FileName_DAPI]]
    info = [split(ii, "_") for ii in Array(ima[!, :Image])]    
    ima[!, :Timepoint] = extract_metainformation(info, 
        ["V 0", "V 1"], 
        [0, 120])  
    ima[!, :Plate] = [ii[4] for ii in info]

    # Describe samples present in control plate
    ima_controls = ima[ima[!, :Plate].== "Controls", :]
    ima_controls2 = ima[ima[!, :Plate].== "Control", :]
    ima_controls = join_in_all_common_columns(ima_controls, ima_controls2)
    info_controls  = [split(ii, "_") for ii in Array(ima_controls[!, :Image])]   
    ima_controls[!,:Genotype] = extract_metainformation(info_controls, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Ifnb1KO", "WT","Rad21KO", "Ifnb1KO"])
    ima_controls[!,:Probe] = extract_metainformation(info_controls, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Cxcl10", "Cxcl10","Cxcl10", "Mx1"])
    
     # Describe samples present in mix plate
    ima_mix = ima[ima[!,:Plate].== "Mix", :]
    info_mix  = [split(ii, "_") for ii in Array(ima_mix[!, :Image])]   
    ima_mix[!,:Genotype] = extract_metainformation(info_mix, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Ifnb1KO+WT", "Ifnb1KO+Rad21KO","Ifnb1KO+WT", "Ifnb1KO+Rad21KO"])
    ima_mix[!, :Probe] = extract_metainformation(info_mix, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Cxcl10", "Cxcl10","Mx1", "Mx1"])
        
    imdata = join_in_all_common_columns(ima_mix, ima_controls)
    
     imdata[!,:Sample] = [imdata[!, :Genotype][ii]*"_"*string(imdata[!, :Timepoint][ii])*"_"*string(imdata[!, :Probe][ii]) for ii in 1:nrow(imdata)]

    
    csv = CSV.read(cell_file, DataFrames.DataFrame)
    csv[!, :Rep] = [1 for i in 1:nrow(csv)]
    csv = innerjoin(csv, imdata, on = :ImageNumber)
    
    csv[!, :Image_Cell] = [csv[ii, :Image]*"__Cell_CP_"*string(csv[ii, :ObjectNumber]) for ii in 1:nrow(csv)]
    return csv
end


function CellInfo4(dir)
    image_file = get_files_ending_with(dir, "Image.csv")
    cell_file = get_files_ending_with(dir, "TrueCells.csv")

    if .&(length(cell_file) == 1,  length(image_file) == 1) # Make sure that there is only one file with the image information
        image_file = string(dir, image_file[1])
        cell_file = string(dir, cell_file[1])
    else
        error("Image.csv or Cell file not found")
    end
    
     ima = CSV.read(image_file, DataFrames.DataFrame)[!, [:FileName_DAPI, :ImageNumber]]
     ima[!, :Image] = [split(split(ii, "_C1")[1], "_MAX")[1] for ii in ima[!, :FileName_DAPI]]
    info = [split(ii, "_") for ii in Array(ima[!, :Image])]    
    ima[!, :Timepoint] = extract_metainformation(info, 
        ["V 0", "V 1"], 
        [0, 120])  
    ima[!, :Plate] = [ii[3] for ii in info]

    # Describe samples present in control plate
    ima_controls = ima[ima[!, :Plate].== "Controls", :]
    info_controls  = [split(ii, "_") for ii in Array(ima_controls[!, :Image])]   
    ima_controls[!, :Genotype] = extract_metainformation(info_controls, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Ifnb1KO", "Rad21KO","WT", "Ifnb1KO"])
    ima_controls[!, :Probe] = extract_metainformation(info_controls, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Cxcl10", "Cxcl10","Cxcl10", "Mx1"])
    
     # Describe samples present in mix plate
    ima_mix = ima[ima[!, :Plate].== "Mix", :]
    info_mix  = [split(ii, "_") for ii in Array(ima_mix[!, :Image])]   
    ima_mix[!, :Genotype] = extract_metainformation(info_mix, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Ifnb1KO+WT", "Ifnb1KO+Rad21KO","Ifnb1KO+WT", "Ifnb1KO+Rad21KO"])
    ima_mix[!, :Probe] = extract_metainformation(info_mix, 
        ["U 0", "U 1", "U 2", "U 3"], 
        ["Cxcl10", "Cxcl10","Mx1", "Mx1"])
        
    imdata = join_in_all_common_columns(ima_mix, ima_controls)
    
     imdata[!, :Sample] = [imdata[!,:Genotype][ii]*"_"*string(imdata[!,:Timepoint][ii])*"_"*string(imdata[!,:Probe][ii]) for ii in 1:nrow(imdata)]

    
    csv = CSV.read(cell_file, DataFrame)
    csv[!, :Rep] = [3 for i in 1:nrow(csv)]
    csv = innerjoin(csv, imdata, on = :ImageNumber)
    
    csv[!, :Image_Cell] = [csv[ii, :Image]*"__Cell_CP_"*string(csv[ii, :ObjectNumber]) for ii in 1:nrow(csv)]
    return csv
end




function FQ_summary_MATURE(dir)
    file = get_files_containing(dir, "summary_MATURE")
    n = length(file)
    if n != 1
        error("ERROR: $n files found")
    end
    filename = string(dir, file[1])
    
     df = CSV.read(filename,DataFrame, delim = '\t', header = 5, skipto = 6)
    rename!(df, :CELL => :Cell)
    rename!(df, :FILE => :Image)
    
    df = FQfiles.fix_image_name(df)
    df = column_fusion(df, :Image, :Cell)

end


root1 = "/Volumes/lymphdev\$/IreneR/Confocal/CellMix_Mx1_Cxcl10/Rep1/"
root2 = "/Volumes/lymphdev\$/IreneR/Confocal/CellMix_Mx1_Cxcl10/Rep2/"
root3 = "/Volumes/lymphdev\$/IreneR/Confocal/CellMix_Mx1_Cxcl10/Rep3/"
root4 = "/Volumes/lymphdev\$/IreneR/Confocal/CellMix_Mx1_Cxcl10/Rep4/"

gene = "Cxcl10"
tss_c4 = root1 * "Segmentation_C2/_FQ_outline/_TS_detect"
imagesfolder = root1 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp1.csv", d)

tss_c4 = root2 * "Segmentation_C2/_FQ_outline/_TS_detect"
imagesfolder = root2 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp2.csv", d)

tss_c4 = root3 * "Segmentation_C2/_FQ_outline/_TS_detect"
imagesfolder = root3 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp3.csv", d)

tss_c4 = root4 * "Segmentation_C2/_FQ_outline/_TS_detect"
imagesfolder = root4 * "/tiff3D"
images_pat = TSSs.get_image_patterns(tss_c4)
d = TSSs.TSS_raw_quant(tss_c4, tss_c4, imagesfolder, 2; xy = 0.189, zx = 0.5)
CSV.write("TSS_raw/"*gene*"_exp4.csv", d)


gene = "Cxcl10"

tss1 = CSV.read("TSS_raw/"*gene*"_exp1.csv", DataFrames.DataFrame)
tss2 = CSV.read("TSS_raw/"*gene*"_exp2.csv", DataFrames.DataFrame)
tss3 = CSV.read("TSS_raw/"*gene*"_exp3.csv", DataFrames.DataFrame)
tss4 = CSV.read("TSS_raw/"*gene*"_exp4.csv", DataFrames.DataFrame)

dot1_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 1)
dot2_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 1)
dot3_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 1)
dot4_r1 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root4 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 1)

dot1_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 2)
dot2_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 2)
dot3_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 2)
dot4_r2 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root4 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 2)

dot1_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root1 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 3)
dot2_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root2 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 3)
dot3_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root3 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 3)
dot4_r3 = TSSs.int_brightest_pixel( TSSs.read_tiff_as_gray(root4 * "/type6_TSS/_mRNA_AVG_ns.tif").*nor; radious = 3)

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



CSV.write("TSS_avgdot/"*gene*"_exp1.csv", tss1)
CSV.write("TSS_avgdot/"*gene*"_exp2.csv", tss2)
CSV.write("TSS_avgdot/"*gene*"_exp3.csv", tss3)
CSV.write("TSS_avgdot/"*gene*"_exp4.csv", tss4)

cp_dir = normpath(ENV["Code"], "Databases2/smFISH_Data/CellMix1_Mx1_Cxcl10_1/CP_results/")
cells1 = CellInfo1(cp_dir);
cp_dir = normpath(ENV["Code"], "Databases2/smFISH_Data/CellMix1_Mx1_Cxcl10_2/CP_results/")
cells2 = CellInfo2(cp_dir);
cp_dir = normpath(ENV["Code"], "Databases2/smFISH_Data/CellMix1_Mx1_Cxcl10_3/CP_results/")
cells3 = CellInfo(cp_dir);
cp_dir = normpath(ENV["Code"], "Databases2/smFISH_Data/CellMix1_Mx1_Cxcl10_4/CP_results/")
cells4 = CellInfo4(cp_dir);

exp1 = CSV.read("TSS_avgdot/Cxcl10_exp1.csv", DataFrames.DataFrame)
exp1[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in exp1[!,:Image]]
CSV.write("Cxcl10/Cxcl10_exp1_TSS.csv", exp1)
CSV.write("Cxcl10/Cxcl10_exp1_CP.csv", cells1)

exp1 = CSV.read("TSS_avgdot/Cxcl10_exp2.csv", DataFrames.DataFrame)
exp1[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in exp1[!,:Image]]

CSV.write("Cxcl10/Cxcl10_exp2_TSS.csv", exp1)
CSV.write("Cxcl10/Cxcl10_exp2_CP.csv", cells2)

exp1 = CSV.read("TSS_avgdot/Cxcl10_exp3.csv", DataFrames.DataFrame)
exp1[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in exp1[!,:Image]]

CSV.write("Cxcl10/Cxcl10_exp3_TSS.csv", exp1)
CSV.write("Cxcl10/Cxcl10_exp3_CP.csv", cells3)

exp1 = CSV.read("TSS_avgdot/Cxcl10_exp4.csv", DataFrames.DataFrame)
exp1[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in exp1[!,:Image]]

CSV.write("Cxcl10/Cxcl10_exp4_TSS.csv", exp1)
CSV.write("Cxcl10/Cxcl10_exp4_CP.csv", cells4)

fq_dir1 = "/Volumes/lymphdev\$/IreneR/Confocal/CellMix_Mx1_Cxcl10/Rep1/type6_TSS/"
fq_dir2 = "/Volumes/lymphdev\$/IreneR/Confocal/CellMix_Mx1_Cxcl10/Rep2/type6_TSS/"
fq_dir3 = "/Volumes/lymphdev\$/IreneR/Confocal/CellMix_Mx1_Cxcl10/Rep3/type6_TSS/"
fq_dir4 = "/Volumes/lymphdev\$/IreneR/Confocal/CellMix_Mx1_Cxcl10/Rep4/type6_TSS/"



fq1 = FQ_summary_MATURE(fq_dir1)
fq2 = FQ_summary_MATURE(fq_dir2)
fq3 = FQ_summary_MATURE(fq_dir3)
fq4 = FQ_summary_MATURE(fq_dir4)


CSV.write("Cxcl10/Cxcl10_exp1_FQ.csv", fq1)
CSV.write("Cxcl10/Cxcl10_exp2_FQ.csv", fq2)
CSV.write("Cxcl10/Cxcl10_exp3_FQ.csv", fq3)
CSV.write("Cxcl10/Cxcl10_exp4_FQ.csv", fq4)



