mismatchPlotRegion <- function(files, region, win = 25, ref, ...) {

    if(width(region) != 1)
        warning("Width of region must be 1.")
    
    n = length(files)
    
    ## pseudo sample data
    if(is.null(names(files))) {
        names = sub(".bam$", "", basename(files))
    }
    sd = data.frame(Type = factor(letters[1:n]), Column = 1:n, Sample = factor(names), Patient = factor("1"))
    
    region = resize(region, width = 2*win+1, fix = "center")
    ref_seq = suppressWarnings(getSeq(ref, region))[[1]]
    theData = lapply(files, function(bamf){
        tallyBAM(file = bamf, chr = as.character(seqnames(region)), start = start(region), stop = end(region), reference = ref_seq, ...)
    })
    theMergedData = mergeTallies(theData)
    
    possible_error = tryCatch(p <- mismatchPlot( theMergedData, sd, samples=sd$Sample, win, win+1, plotReference = FALSE ), error = function(e) {e})

    if(inherits(possible_error, "error"))
        return(NULL)
    else
        return(p)
}


tallyPlot <- function(file, region, ref, nCycles = 0, minQual = 0, minFreq = 0, ...) {

    n = width(region)
    
    tally = tallyBamRegion(file, region, nCycles, minQual)
    tally = selectStrand(tally, "both")
    colnames(tally) = c("A", "C", "G", "T")
    
    ref = ds2int(getSeq(ref, region)[[1]])
    
    ## add 'N/R' base to tally
    tally = cbind(tally, N = tally[cbind(1:n, ref)])
    tally[cbind(1:n, ref)] = 0
    
    d = melt(tally)
    #dn = subset(d, value > 0) ## causes NOTE in check
    dn = d[d$value > 0, ]
    
    gn = GRanges(seqnames(region),
        IRanges(start(region):end(region), width = 1))
    gx = gn[dn$Var1]
    gx$base = dn$Var2
    gx$count = dn$value
    #gx = sort(gx)
        
    p = autoplot(gx,
        aes_string(x = "start", y = "count", fill = "base"),
        stat = "identity",
        position = "stack",
        geom = "bar",
        width = 1, ...)
    p = p + baseFillScale() + theme_bw()
    
    return(p)
}
