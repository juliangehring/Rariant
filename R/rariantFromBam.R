rariantFromBam <- function(test, control, region, beta = 0.95, alpha = 1 - beta, select = TRUE, consensus, resultFile, strand = c("both", "plus", "minus"), nCycles = 10, minQual = 20, block = 1e4, value = TRUE, criteria = c("both", "any", "fet", "ci"), verbose = 1e7) {

    ## input arguments
    args = as.list(match.call())[-c(1:4)] ## ignore the function name and 'test', 'control', 'region'
    strand = match.arg(strand)
    criteria = match.arg(criteria)

    write_file = !missing(resultFile)

    ## test the input files
    seq_info_test = seqinfo(test)
    seq_info_control = seqinfo(control)
    if(!identical(seq_info_test, seq_info_control))
        warning("The BAM headers of 'test' and 'control' differ!")

    suppressWarnings(seqinfo(region) <- seq_info_control[seqlevels(region)])
    region = trim(region)
    chunks = region
    if(max(width(chunks)) > block) ## if not, result is the same, but we save time
        chunks = subdivideGRanges(chunks, block) ## tested

    #w = width(chunks)
    #d = diff(w)
    
    header = TRUE
    val = list()
    for(i in seq_along(chunks)) {
        ## the current region
        roi = chunks[i]
        n = width(roi)

        test_tally = tallyBamRegion(test, roi, nCycles, minQual)
        control_tally = tallyBamRegion(control, roi, nCycles, minQual)

        rar = rariantFromMatrix(test_tally, control_tally, roi, beta, alpha, select, consensus, strand, criteria)

        ## write output file
        if(write_file && nrow(rar) > 0) {
            write.table(rar, resultFile, append = !header, quote = FALSE, sep = "\t", row.names = FALSE, col.names = header)
            header = FALSE
        }
        
        if(value) {
            val[[i]] = rar
        }     
    }
    
    res = NULL ## must be set
    if(value & length(val) > 0) {
        res = rbind_all(val) ## or unlist on a GRL
        #if(!is.null(res)) { ## TODO: needed?
        res = df2gr(res)
        #}
    } else {
        res = GRanges()
    }
    ## Setting the seqinfo here means that there is none
    ## if no variant is returned
    seqinfo(res) = seq_info_test[seqlevels(res)]
    metadata(res) = args
    
    return(res)
}
