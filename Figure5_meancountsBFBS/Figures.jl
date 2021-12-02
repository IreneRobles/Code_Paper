using PrettyPlotting
function dysplay_correlation2(nonantb, x, y; style = "RNAseqcomp", hue = "GeneSymbol", palette = nothing, xy = (0, 4))

    x1 = [ii for ii in nonantb[!,x]]
    boolx = .! isnan.(x1) .& .! isinf.(x1)
    y1 = [ii for ii in nonantb[!,y]]
    booly = .! isnan.(y1) .& .! isinf.(y1)
    bool = boolx.*booly
    x1 = x1[bool]; y1 = y1[bool]
    
    
    pdt = Pandas.DataFrame(sort(nonantb[bool, :], :Genotype, rev = true))

    
    py"""
    import seaborn as sns
    sns.lmplot(data = $pdt, x= $x, y = $y, hue = $hue, fit_reg = 1, height=3,aspect=1, palette = $palette)
    """
    
    pretty_axes2()
    test = R"cor.test($x1,$y1)"
    pval = round(test[3][1], sigdigits = 2)
    corr = round(test[4][1], sigdigits = 2)
    annotate("""
        r = $corr
        p-value = $pval
        """, xy = xy, va = "center"
    )  
    
    end 