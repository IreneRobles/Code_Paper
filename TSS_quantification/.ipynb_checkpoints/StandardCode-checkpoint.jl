ENV["Code"] = "../../Code"
for folder in readdir(ENV["Code"]); push!(LOAD_PATH, normpath(ENV["Code"], folder));end
include("../../Code/TSSs/src/TSSs.jl")
using ProgressMeter
using DataFrames
using Distances
using CSV
using NoLongerProblems_FileHandling
using NoLongerProblems
using NoLongerProblems_Pandas
using DataFrames
using FQfiles
using HypothesisTests
using MultipleTesting
using Seaborn
import Pandas

using PyPlot
using PrettyPlotting


function add_probetypes(e)
    col = Symbol("Probe (type)")
    type1 = []; type4 = []; type6 = []; 
    
    for ii in e[!,col]
        probes = split(replace(ii, " " => ""), ",")
        probes_pretty = [split(jj, "(")[1] for jj in probes]
        
        t1_ind = findall(x -> occursin("(1)", x), probes)
        t4_ind = findall(x -> occursin("(4)", x), probes) 
        t6_ind = findall(x -> occursin("(6)", x), probes) 
        
        t1 = if length(t1_ind) == 1 probes_pretty[t1_ind][1] else "NA" end; 
        t4 = if length(t4_ind) == 1 probes_pretty[t4_ind][1] else "NA" end; 
        t6 = if length(t6_ind) == 1 probes_pretty[t6_ind][1] else "NA" end; 
        
        push!(type1, t1); push!(type4, t4); push!(type6, t6); 
    end
    e[!,:type1] = type1; e[!,:type4] = type4; e[!,:type6] = type6; 
    e
    
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
    
    ima = DataFrames.DataFrame(CSV.read(image_file, DataFrame))[:, [:FileName_DAPI, :ImageNumber]]
    ima[!,:Image] = [split(split(ii, "_C1")[1], "_MAX")[1] for ii in ima[!,:FileName_DAPI]]
    
    ima[!,:WELL] = [split(ii, " ") for ii in ima[!,:Image]]
    ima[!,:WELL] = [1 for ii in ima[!,:WELL]]
    
     ima[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in ima[!,:FileName_DAPI]]
    
    exp_df[!,:Well] = [split(ii, " (")[1] for ii in exp_df[!,:Well]]

    im_ = innerjoin(exp_df, ima, on = :Well, makeunique=true)
        

    
   
     cells = DataFrames.DataFrame(CSV.read(cell_file, DataFrame))
    cells = innerjoin(im_, cells, on = :ImageNumber)
    
    cells[!,:Image_Cell] = [cells[ii, :Image]*"__Cell_CP_"*string(cells[ii, :ObjectNumber]) for ii in 1:nrow(cells)]
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
    
    ima = DataFrames.DataFrame(CSV.read(image_file, DataFrame))[:, [:FileName_DAPI, :ImageNumber]]
    ima[!,:Image] = [split(split(ii, "_C1")[1], "_MAX")[1] for ii in ima[!,:FileName_DAPI]]
    
    ima[!,:WELL] = [split(ii, " ") for ii in ima[!,:Image]]
    ima[!,:WELL] = [1 for ii in ima[!,:WELL]]
    
     ima[!,:Well] = [split(split(ii, "S 0_")[2], "_X")[1] for ii in ima[!,:FileName_DAPI]]
    
    exp_df[!,:Well] = [split(ii, " (")[1] for ii in exp_df[!,:Well]]

    im_ = innerjoin(exp_df, ima, on = :Well, makeunique=true)
        

    
   
     cells = DataFrames.DataFrame(CSV.read(cell_file))
    cells = innerjoin(im_, cells, on = :ImageNumber)
    
    cells[:Image_Cell] = [cells[ii, :Image]*"__Cell_CP_"*string(cells[ii, :ObjectNumber]) for ii in 1:nrow(cells)]
    return dropmissing(cells)
    
end



function FQ_summary_MATURE(dir)
    file = get_files_containing(dir, "summary_MATURE")
    n = length(file)
    if n != 1
        error("ERROR: $n files found")
    end
    filename = string(dir, file[1])
    
     df = CSV.read(filename, DataFrame, delim = '\t', header = 5, skipto = 6)
    rename!(df, :CELL => :Cell)
    rename!(df, :FILE => :Image)
    
    df = FQfiles.fix_image_name(df)
    df = column_fusion(df, :Image, :Cell)

end

function cells_per_sample(df)
    sams = unique(df[!,:Name])
    ncells = [ count(x -> x == sam, df[!,:Name]) for sam in sams]
    new_df = DataFrames.DataFrame(:Sample=> sams, :N_Cells=> ncells)
end


function assign_red_green_labels(REP, max_green_for_red, min_red_for_red, max_red_for_green,  min_green_for_green, nrep)
    function DrawRectangles()
        topcorner = 1
        plot([max_green_for_red, max_green_for_red], [min_red_for_red, topcorner], c = "red")
        plot([0, 0], [min_red_for_red, topcorner], c = "red")
        plot([0, max_green_for_red], [min_red_for_red, min_red_for_red], c = "red")
        plot([0, max_green_for_red], [1, 1], c = "red")

        plot([min_green_for_green, 1], [max_red_for_green, max_red_for_green], c = "green")

        plot([1, 1], [0, max_red_for_green], c = "green")

        plot([min_green_for_green, min_green_for_green], [max_red_for_green, 0], c = "green")
        plot([min_green_for_green, 1], [0, 0], c = "green")
    end
    
    
    
    
    green = REP[:Intensity_UpperQuartileIntensity_FilteredGreen]
    red = REP[:Intensity_UpperQuartileIntensity_FilteredRed]
    labels = ["NotKnown" for ii in 1:length(red)]

    for ii in 1:length(red)
        r = red[ii]; g = green[ii]

        if .&(g > min_green_for_green, r < max_red_for_green)
            labels[ii] = "Green"
        elseif .&(r > min_red_for_red, g < max_green_for_red)
            labels[ii] = "Red"
        end

    end

    REP[:Label] = labels

    redcount = count(x -> x == "Red", labels)
    println("Red = $redcount")

    gcount = count(x -> x == "Green", labels)
    println("Green = $gcount")

    REP[:UpperQuartileIntensity_Green] =  REP[:Intensity_UpperQuartileIntensity_FilteredGreen] 
    REP[:UpperQuartileIntensity_Red] =  REP[:Intensity_UpperQuartileIntensity_FilteredRed]

    y = :UpperQuartileIntensity_Red
    x = :UpperQuartileIntensity_Green
    figure(figsize = (10, 5))

    subplot(1, 2, 1)

    sam = "WT"
    df = REP[REP[:Cells].== sam, :]
    scatter(df[x], df[y], s = 1, label = sam, c = "blue")
    
     sam = "Rad21KO"
    df = REP[REP[:Cells].== sam, :]
    scatter(df[x], df[y], s = 1, label = sam, c = "orange")
    
    
    
    DrawRectangles(); pretty_axes2(); ylabel(y); xlabel(x)
    
    subplot(1, 2, 2)

    sam = "Rad21KO+WT"
    df = REP[REP[:Cells].== sam, :]
    scatter(df[x], df[y], s = 1, label = sam, c = "purple")
    
    
    DrawRectangles(); pretty_axes2(); ylabel(y); xlabel(x)
    
    return REP

end

function TSSs.TSS_raw_quant(t2, tss_folder, image_folder, n; xy = 0.189, zx = 0.5)
    
    images_pat = TSSs.get_image_patterns(t2)
    
    p = Progress(length(images_pat), 1)

    dfs = []

    for a in 1:length(images_pat)
        next!(p)
        

        pat = images_pat[a]

        # Get all the TSS in the images
        a = TSSs.find_outline(tss_folder, pat)
        tss = TSSs.cell_tss_dict(a)
        tss = TSSs.meassure_tss(tss, image_folder, pat, n, xy = xy, zx = zx)
    
        df = TSSs.analysis_singlechannel(tss; image = pat,  tss_name = :TSS2)
        push!(dfs, df)

    end
    
    d = TSSs.join_in_all_common_columns(dfs...)
    return d
    
end

function bs_bf(exp; limit = 0.1, r = "r1")
    sams = unique(exp[!,:Name])
    df_sp = split_by(exp, :Name)
    ncells = [ count(x -> x == sam, exp[!,:Name]) for sam in sams]
    #bursts1 = [ count(x -> x > limit, df_sp[sam][!,Symbol(string("TSS1_", r))]) for sam in sams]
    #bursts2 = [ count(x -> x > limit, df_sp[sam][!,Symbol(string("TSS2_", r))]) for sam in sams]
    bursts1 = [ count(x -> x > limit, df_sp[sam][!,Symbol("locus1_int2_TSS2")]) for sam in sams]
    bursts2 = [ count(x -> x > limit, df_sp[sam][!,Symbol(string("locus2_int2_TSS2"))]) for sam in sams]
    new_df = DataFrames.DataFrame(:Sample=> sams, :N_Cells=> ncells, :N_TSS1 => bursts1, :N_TSS2 => bursts2)
    new_df[!,:BF] = (new_df[!,:N_TSS1] .+ new_df[!,:N_TSS2])./2new_df[!,:N_Cells]
    new_df
end
