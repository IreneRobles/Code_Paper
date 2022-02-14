using HypothesisTests

function anova(a, b)
    taga = ["A" for ii in 1:length(a)]
    tagb = ["B" for ii in 1:length(b)]
    my_data = DataFrames.DataFrame(
        weight = append!(a, b),
        group = append!(taga, tagb)
    )
    
    R"""
    # Compute the analysis of variance
    res.aov <- aov(weight ~ group, data = $my_data)
    # Summary of the analysis
    summary(res.aov)[[1]][["Pr(>F)"]][1]
    """[1]
end

function anova(a, b, repa, repb)
    taga = ["A" for ii in 1:length(a)]
    tagb = ["B" for ii in 1:length(b)]
    my_data = DataFrames.DataFrame(
        weight = append!(a, b),
        group = append!(taga, tagb),
        rep = append!(repa, repb)
    )
    
    R"""
    # Compute the analysis of variance
    res.aov <- aov(weight ~ group , data = $my_data)
    # Summary of the analysis
    summary(res.aov)[[1]][["Pr(>F)"]][1]
    """[1]
end

function do_mantelhaen(df, comp1, comp2)
    s1 = df[df[!,:Sample] .== comp1, :]
    s2 = df[df[!,:Sample] .== comp2, :]
    tb_test = join_in_all_common_columns(s1, s2)
    tb1 = tb_test; tb1[!,:Burst] = ["Yes" for ii in 1:nrow(tb_test)]
    tb1[!,:Count] = tb1[!,:nburst]
    tb1 = tb1[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb2 = tb_test; tb2[!,:Burst] = ["No" for ii in 1:nrow(tb_test)]
    tb2[!,:Count] = tb2[!,:n_cells] * 2 - tb2[!,:nburst]
    tb2 = tb2[!, [:Sample, :Count, :Rep, :Burst]]
    
    tb = join_in_all_common_columns(tb1, tb2)

    
test = R"""
#library("psych")
#library("vcd")
#library("DescTools")
#library("rcompanion")
library("dplyr")
options(warn=-1)

tb = $tb
tb$Count = as.integer(tb$Count)
    
Data <- mutate(tb,
           Sample = factor(Sample, levels=unique(Sample)),
           Burst = factor(Burst, levels=unique(Burst)),
           Rep = factor(Rep, levels=unique(Rep))
           )

# Last variable is the strata (the variable that is not check for assotiation)
Data.xtabs <- xtabs(Count ~ Burst + Sample + Rep, 
                       data=Data)

ftable(Data.xtabs)  

mantelhaen.test(Data.xtabs)
"""

return test[3][1]
end

