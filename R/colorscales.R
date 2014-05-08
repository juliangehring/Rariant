eventFillScale <- function() {
    values = c(somatic = "#fc8d62", hetero = "#e78ac3", undecided = "#8da0cb", powerless = "#a6d854", none = "white")
    return(scale_fill_manual(values = values))
}

verdictColorScale <- function() {
    values = c(present = "#fc8d62", inbetween = "#8da0cb", dontknow = "#a6d854", absent = "lightgray")
    return(scale_color_manual(values = values))
}

baseFillScale <- function() {
    ##library(RColorBrewer)
    ##cols = brewer.pal(8, "Set2")[c(1, 2, 3, 6)] ## 8 for gray
    ##cols = c(cols, "lightgray")
    ##names(cols) = c("A", "C", "G", "T", "R")
    values = c(A = "#57D283", C = "#FC8D62", G = "#5F82D1", T = "#FFD92F", N = "lightgray")
    return(scale_fill_manual(values = values))
}

rateFillScale <- function() {
    return(scale_fill_gradient2(limits = c(-1, 1)))
}
