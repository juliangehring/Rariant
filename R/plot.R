plotConfidenceIntervals <- function (x, ylim = c(-1.05, 1.05), color = NULL, ...) {

    dfh = data.frame(h = -1:1)
    
    ## convert until we have ggbio support for VRanges
    if(inherits(x, "GRanges"))
        x = vr2gr(x)

    if ( inherits(x, "GRanges") ) {
        p = autoplot(x, geom = "pointrange", aes_string(y = "d", ymin = "lower", ymax = "upper", color = color), ...)
    }

    if ( inherits(x, "data.frame") ) {
        if ( is.null(x$start) )
            x$start = 1:nrow(x)

        p = ggplot(x) + geom_pointrange(aes_string(x = "start", y = "d", ymin = "lower", ymax = "upper", color = color), ...)
    }

    p = p + geom_hline(aes_string(yintercept = "h"), data = dfh, color = "darkgray", linetype = "dashed") + theme_bw() + coord_cartesian(ylim = ylim) + xlab("Genomic position") + ylab("Variant Rates Difference")
        
    return(p)
}


plotAbundanceShift <- function (x, ylim = c(-0.05, 1.05), ...) {

    dfh = data.frame(h = 0:1)
  
    ## convert until we have ggbio support for VRanges
    if(inherits(x, "GRanges"))
        x = vr2gr(x)

    x$bottom = pmin(x$p1, x$p2)
    x$top = pmax(x$p1, x$p2)
    x$shift = factor(ifelse(x$p1 < x$p2, "loss", "gain"), levels = c("gain", "loss"))
    
    if ( inherits(x, "GRanges") ) {
        p = autoplot(x, geom = "linerange", aes_string(ymin = "bottom", ymax = "top", color = "shift"), ...)
    }
    
    if ( inherits(x, "data.frame") ) {
        if ( is.null(x$start) )
            x$start = 1:nrow(x)
        
        p = ggplot(x) + geom_linerange(aes_string(x = "start", ymin = "bottom", ymax = "top", color = "shift"), ...) + geom_point(aes_string(x = "start", y = "p1"), ...) + geom_point(aes_string(x = "start", y = "p2"), ...)
    }
    
    p = p + geom_hline(aes_string(yintercept = "h"), data = dfh, color = "darkgray", linetype = "dashed") + theme_bw() + coord_cartesian(ylim = ylim) + xlab("Genomic position") + ylab("Non-consensus rates")
    
    return(p)
}


plotShiftSupport <- function (x, group = NULL, alpha = 0.5, size = 2)  {

    df = data.frame(vf = x$d, cov = x$testDepth)
    if (!is.null(group)) 
        df[, group] = as(x, "data.frame")[, group]
    df = df[!is.na(df$vf), ]
    p = ggplot(df) + geom_point(aes_string(x = "vf", y = "cov", 
        col = group), size = size, alpha = alpha) + xlim(0, 1) + 
        xlab("Shift") + ylab("Test sequencing depth") + theme_bw() + 
        scale_y_log10()
    
    return(p)
}
