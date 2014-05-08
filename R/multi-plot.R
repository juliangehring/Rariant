evidenceHeatmap <- function(x, fill = "d", color = "outside",
                           height = 0.9, width = 0.95, size = 1.5,
                           xvar = "sample", yvar = "loc", ...) {

    mc = mergeCalls(x)

    p = ggplot(mc) +
        geom_tile(aes_string(x = xvar, y = yvar, fill = fill, color = color),
                  size = size, height = height, width = width, ...) +
                      theme_bw() + xlab("Sample") + ylab("Variant")
    
    return(p)
}
