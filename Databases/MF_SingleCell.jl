using DataFrames #RNAseq_v2
import Distances
using NoLongerProblems_FileHandling
using NoLongerProblems
using CSV

export singlecell_distancearray, get_cell_to_cell_distances_as_vector, get_cell_to_cell_distances
export SergiSingleCell_colData, SergiSingleCell_counts
export get_sample, get_Bhatt_genes

function get_sample(sampl; counts = cd, colData = coldata)
    cols = [Symbol(i)  for i in collect(coldata[coldata[!,:Sample].==sampl, :RowName])]
    counts[:, append!(cols, [:GeneSymbol])]
end

cd = DataFrame!(CSV.File(string("../Databases/MF_SingleCell/data_for_julia/","cd.csv")));
coldata = DataFrame!(CSV.File(string("../Databases/MF_SingleCell/data_for_julia/","coldata.csv")));
coldata[!, :Timepoint] = [split(i, "_")[end] for i in coldata[!,:Sample]]
coldata

df_summary = DataFrames.DataFrame()
df_summary[!,:Sample] = unique(coldata[!,:Sample])
df_summary[!,:Timepoint] = [split(i, "_")[end] for i in df_summary[!,:Sample]]
df_summary[!,:Genotype] = [split(i, "_")[1] for i in df_summary[!,:Sample]]
df_summary[!,:Cells] = [size(get_sample(string(i)))[2] for i in df_summary[!,:Sample]]

df_summary

function SergiSingleCell_colData()
    return coldata
end

function SergiSingleCell_summary()
    return df_summary
end


function SergiSingleCell_counts()
    return cd
end

function get_Bhatt_genes(; counts = cd, class = unique(Bhatt2012.inducible_genes_figure3()[:Class]))
    bhatt = Bhatt2012.inducible_genes_figure3()
    f(x) = in(x, class)
    bhattgenes = bhatt[[f(i) for i in bhatt[!,:Class]], :]
    selected_genes = bhattgenes[!, [:GeneSymbol]]
    return join(counts, selected_genes, on = :GeneSymbol)
end

function get_Bhatt_LPS_Induced_genes(; counts = SergiSingleCell.cd, group = 1:3)
    bhatt = Bhatt2012.all_genes_figure2()
    f(x) = in(x, group)
    bhattgenes = bhatt[[f(i) for i in bhatt[!,:Group]], :]
    selected_genes = bhattgenes[:, [:GeneSymbol]]
    return join(counts, selected_genes, on = :GeneSymbol)
end

function get_Bhatt_LPS_Silenced_genes(; counts = SergiSingleCell.cd, group = 4:6)
    bhatt = Bhatt2012.all_genes_figure2()
    f(x) = in(x, group)
    bhattgenes = bhatt[[f(i) for i in bhatt[!,:Group]], :]
    selected_genes = bhattgenes[:, [:GeneSymbol]]
    return join(counts, selected_genes, on = :GeneSymbol)
end

function get_Bhatt_LPS_Response_genes(; counts = SergiSingleCell.cd, group = 1:6)
    bhatt = Bhatt2012.all_genes_figure2()
    f(x) = in(x, group)
    bhattgenes = bhatt[[f(i) for i in bhatt[!,:Group]], :]
    selected_genes = bhattgenes[:, [:GeneSymbol]]
    return join(counts, selected_genes, on = :GeneSymbol)
end


function get_RNAseq_DE_genes(; counts = cd, padj = 0.05, limit_expression = 100)
    rnaseq = SergiData.summary_SergiMF_data()
    gene_selection_rnaseq = RNAseq_v2.genes_DE_in_one_condition_at_least(rnaseq, limit = padj)
    gene_selection_rnaseq = RNAseq_v2.genes_expressed_in_one_condition_at_least(gene_selection_rnaseq, limit = limit_expression)
    selected_genes = gene_selection_rnaseq[!, [:GeneSymbol]]
    return join(counts, selected_genes, on = :GeneSymbol)
end

function get_RNAseq_DE_genes_in_condition(condition_padj; counts = SergiSingleCell.cd, padj = 0.05)
    rnaseq = SergiData.summary_SergiMF_data()
    gene_selection_rnaseq = rnaseq[rnaseq[condition_padj] .< padj, :]
    selected_genes = gene_selection_rnaseq[:, [:GeneSymbol]]
    return join(counts, selected_genes, on = :GeneSymbol)
end





function singlecell_distancearray(df; distance = Distances.euclidean)
    d_1 = df[setdiff(names(df), [:GeneSymbol])]
    cells = ncol(d_1)
    cell_cell_distances = Array{Any, 2}(cells, cells)
    for i in 1:cells
        cell1 = d_1[:, i]
        for j in 1:cells
            cell2 = d_1[:, j]
            cell_cell_distances[i, j] = distance(cell1, cell2)
        end
    end       
    cell_cell_distances
end

function get_cell_to_cell_distances_as_vector(array)
    if size(array)[1] == size(array)[2]
        n_cells = size(array)[1]
        
        vector = Vector{Float64}()
        
        for i in 1:n_cells
            for j in (i+1):n_cells # i + 1 because we do not want to count the distance between a cell and itself
                push!(vector, array[i, j])
            end
        end
        return vector
    else
        "Array must contain the same number of rows and columns"
    end
end

function get_cell_to_cell_distances(df; distance = Distances.euclidean);
    euclideanarray = singlecell_distancearray(df,  distance = distance)
    vector = get_cell_to_cell_distances_as_vector(euclideanarray)
    return vector
end

export get_sample_gene_distribution

function get_sample_gene_distribution(genes; counts = SergiSingleCell_counts(), coldata = SergiSingleCell_colData())
    
    coldata_new_names = [replace(i, ".", "_") for i in coldata[:RowName]]
    coldata[:RowName] = coldata_new_names
    
    for gene in genes
        gene_counts = counts[counts[:GeneSymbol].== gene, :]
        order_list = []

        for col in coldata_new_names
            push!(order_list, gene_counts[1, Symbol(col)])
        end

        coldata[Symbol(gene)] = order_list
    end
    return coldata


end

function single_cell_dot(gene1, gene2)
        
    df = get_sample_gene_distribution([gene1, gene2])
    
    f(x) = in(x, ["WT", "RAD21"])
    
    df = df[[f(i) for i in df[:Genotype]], :]
    
    df = Pandas.DataFrame(df)
    
    lmplot(x=gene1, y=gene2, col="Sample", hue="Genotype", data=df,
           col_wrap=2, palette="muted", size=4)
    
end


function show_expression_replica_genes(genes_...; measure = cor)
    fig = figure(figsize = (8, 2))
    
    genes = [split(i, " ") for i in genes_]
    
    gene_list = Array{String, 1}()
    
    for i in genes
        for j in i
            push!(gene_list, j)
        end
    end
    
    n = length(gene_list)
    
    df = get_sample_gene_distribution(gene_list)
    
    fpkm = SergiData.get_FPKM()
    
    for i in 1:n
        subplot(n, n, i)
        title(gene_list[i])
        SergiData.ReplicateExpression(fpkm, gene_list[i])
        ax = gca()
        if i !=n
            ax[:legend_][:remove]()
        end
        if i == 1
            ylabel("FPKM")
        end
    end
    
    c = 0
    
    plt[:tight_layout]()
    
    correlations_singlecell_genes(genes_..., measure = measure)
end


function correlations_singlecell_genes(genes...; measure = cor)
        
    genes = [split(i, " ") for i in genes]
    
    gene_list = []
    
    for i in genes
        for j in i
            push!(gene_list, j)
        end
    end
    
    n = length(gene_list)
    
    gene_symbols = [Symbol(i) for i in gene_list]
    
    df = DataFrame()
    
    df = get_sample_gene_distribution(gene_list)
    
    f(x) = in(x, ["WT", "RAD21"])
    
    df = df[[f(i) for i in df[:Genotype]], :]
    
    correlations = DataFrame()
    
    correlations[:Sample] = unique(df[:Sample])
    
    
    for i in 1:length(gene_list)
        for j in (i+1):length(gene_list)
            col_ent = Symbol(string("Entropy_", gene_list[i], "_",gene_list[j]))
            col_pear = Symbol(string("Measure_", gene_list[i], "_",gene_list[j]))
            gene2 = gene_symbols[i]
            gene1 = gene_symbols[j]
            
            f(sam, gene) = df[equals(df[:Sample], sam), gene]
            
            cors = [measure(Array{Float64, 1}(f(sam, gene1)), Array{Float64, 1}(f(sam, gene2))) for sam in correlations[:Sample]]
            
            correlations[col_pear] = cors
            
        end
    end
    
    correlations
            
end

function unfold(A)
    V = []
    for x in A
        if x === A
            push!(V, x)
        else
            append!(V, unfold(x))
        end
    end
    V
end

function dropnan(array)
    nan_bool = [isnan(i) == false for i in array]
    array[nan_bool]
end

export calculatefromsinglecelldatatable

function calculatefromsinglecelldatatable(; measure = get_entropy)
     genedata = prepare_genes_for_combined_network(only_protein_coding(one_gene_one_locus_databasetable()))[:, [:EnsemblID, :GeneSymbol]]
    
    
    counts = SergiSingleCell_counts()[:, [:GeneSymbol]]
    
    table = join(genedata, counts, on = :GeneSymbol, kind = :inner)
    
    genes = table[:GeneSymbol]
    
    samples = unique(collect(get_sample_gene_distribution(["Il12b"])[:Sample]))
    
    gene_cors = []
    
    c = 0
 
    for gene in genes
        
        c += 1
        
        if c%500 == 0
            println("$c genes calculated")
        end
        
        g_sym = Symbol(gene) 
        df = get_sample_gene_distribution([gene])
        f(sam, g) = df[equals(df[:Sample], sam), g]
        cors = [measure(Array{Float64, 1}(f(sam, g_sym))) for sam in samples]        
        push!(gene_cors, cors)

    end
    
    for sampl in 1:length(samples)
        sam_sym = Symbol(samples[sampl])
        table[sam_sym] = [i[sampl] for i in gene_cors]
        
    end
    
    return table

end

function calculatefromsinglecelldatatable(; measure = get_entropy, kwargs = Dict())
    
    genedata = prepare_genes_for_combined_network(only_protein_coding(one_gene_one_locus_databasetable()))[:, [:EnsemblID, :GeneSymbol]]
    counts = SergiSingleCell_counts()[:, [:GeneSymbol]]
    table = join(genedata, counts, on = :GeneSymbol, kind = :inner)
    
    genes = table[:GeneSymbol]
    
    samples = unique(collect(get_sample_gene_distribution(["Il12b"])[:Sample]))
    
    gene_cors = []
    
    c = 0
 
    for gene in genes
        c += 1
        if c%500 == 0
            println("$c genes calculated")
        end
        g_sym = Symbol(gene)
        df = get_sample_gene_distribution([gene])
        f(sam, g) = df[equals(df[:Sample], sam), g]
        cors = [measure(Array{Int, 1}(f(sam, g_sym)); kwargs...) for sam in samples]
        push!(gene_cors, cors)
        
    end
    
    for sampl in 1:length(samples)
        sam_sym = Symbol(samples[sampl])
        table[sam_sym] = [i[sampl] for i in gene_cors]
    end
    
    return table
    
end
        
nils_folder = string("../Databases/MF_SingleCell/Nils/AfterRepoolingDKO/")

function read_newNilsDE(filen)
    file = string(nils_folder, filen, ".csv")
    table = CSV.read(file, separator =' ', skipstart = 1, header = false)[:, [3, 4, 5, 6]]
    rename!(table, :x3=> :MeanOveral)
    rename!(table, :x4=> :ResDispDistance)
    rename!(table, :x5=>:ResultDiffResDisp)
    rename!(table, :x6=>:GeneSymbol)
    
    
    group1 = split(filen, "_")[1]
    group2 = split(filen, "_")[2]
    
    table[!,:ResultDiffResDisp] = [replace(i, "Group1", group1) for i in table[!,:ResultDiffResDisp]]
    table[!,:ResultDiffResDisp] = [replace(i, "Group2", group2) for i in table[!,:ResultDiffResDisp]]
    

    result_disp = Symbol(string(:ResultDiffResDisp, "_",filen))
    rename!(table, :ResultDiffResDisp, result_disp)
    
    res_disp = Symbol(string(:ResDispDistance, "_",filen))
    rename!(table, :ResDispDistance, res_disp)
    meanove = Symbol(string(:MeanOveral, "_",filen))
    rename!(table, :MeanOveral, meanove)
    
    table
    
end

function print_total_genes(df)
    n = size(df)[1]
    println("There are $n genes being considered")
end

export genes_more_variable_in
function genes_more_variable_in(df; sampl = "WT", col = :ResultDiffResDisp)
    f(x) = contains(x, sampl)
    df = substitute_NA_for_this(df, col, "NoDiff")
    genes = df[[f(i) for i in df[!,col]], :]
    genes
end


export NillsVariabilityTable_WT_CTCF_DKO
function NillsVariabilityTable_WT_CTCF_DKO()

    # Read files
    ct_dko = read_newNilsDE("CTCF_DKO")
    wt_ct = read_newNilsDE("WT_CTCF")
    wt_dko = read_newNilsDE("WT_DKO")


    # Merge files 

    all = join(ct_dko, wt_ct, on = :GeneSymbol, kind = :outer)
    all = join(all, wt_dko, on = :GeneSymbol, kind = :outer)
    return all
    
end
        

function scde_DE(group1, group2; folder = string("../Databases/MF_SingleCell/scde/DifferentialExpression/LPSinducible"))
    println("Using folder $folder")
    files = readdir_fullpath(folder)
    bool = [contains(i, group1) && contains(i, group2) for i in files]
    files = files[bool]
    if length(files) != 1
        return "Several files found containing both groups"
    end
    
    diffexp = CSV.read(files[1])
    
    rename!(diffexp, :x, :GeneSymbol)
    rename!(diffexp, :mle, :log2FC)
    rename!(diffexp, :lb, :lowerbound)
    rename!(diffexp, :ub, :upperbound)
    
    diffexp[:pvalue] = [z_to_pvalue_2sided(i) for i in diffexp[:Z]]
    diffexp[:padj] = [z_to_pvalue_2sided(i) for i in diffexp[:cZ]]
    
    sort!(diffexp, cols = :log2FC)
    
    return diffexp
 
end

function scde_DE_save(file; folder = folder = string(ENV["Code"], "/Databases/SergiSingleCell/scde/DifferentialExpression/LPSinducible", pvalue = 0.05))
    diffexp = readtable(file)
    
    rename!(diffexp, :x, :GeneSymbol)
    rename!(diffexp, :mle, :log2FC)
    rename!(diffexp, :lb, :lowerbound)
    rename!(diffexp, :ub, :upperbound)
    
    diffexp[:pvalue] = [z_to_pvalue_2sided(i) for i in diffexp[:Z]]
    diffexp[:padj] = [z_to_pvalue_2sided(i) for i in diffexp[:cZ]]
    
    sort!(diffexp, cols = :log2FC)
    

    diffexp = lessthan(diffexp, :padj, pvalue)
    
    writetable(string(folder, "/", split(file, "/")[end]), diffexp)
 
end

function color_samples()
    Dict(
        "WT" => "#4286f4",
        "CTCF" => "#ef4362", 
        "RAD21" => "orange",
        "DKO" => "green",
        "RAD21CTCF_DKO" => "green", 
        "RAD21CTCF" => "green", 
        "UT" => "lightblue",
        "2H" => "blue",
        "8H" => "darkblue",
    )
    
end

