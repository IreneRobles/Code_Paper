using CSV, DataFrames

function bed_writer_cols(tb, filename; chr = :chr, start_ = :start, end_ = :end)
    CSV.write(filename, tb[!,[chr, start_, end_]], delim = '\t', header = false)
end

function bed_readerTADs(file; skipto = 2)
     df = CSV.read(file, DataFrames.DataFrame, header = false, delim = "\t", skipto = skipto)  
    rename!(df, names(df)[1] => :TADchr)
    rename!(df, names(df)[2] => :TADstart)
    rename!(df, names(df)[3] => :TADend)
    rename!(df, names(df)[4] => :TAD)
    df[!,"TAD"] = string.(df[!,"TADchr"], df[!,"TAD"])
    return df
end