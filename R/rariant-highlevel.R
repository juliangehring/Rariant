rariant <- function(test, control, region, beta = 0.95, alpha = 1 - beta, select = TRUE, consensus, resultFile, strand = c("both", "plus", "minus"), nCycles = 10, minMapQual = 20, block = 1e4, value = TRUE, criteria = c("both", "any", "fet", "ci")) {

    ## input arguments
    args = as.list(match.call())[-1] ## ignore the function name
    strand = match.arg(strand)
    criteria = match.arg(criteria)
    #consensus = match.arg(consensus)
    write_file = !missing(resultFile)

    ## test if files exist
    if(!file.exists(test) | !file.exists(control))
        stop("Files 'test' or 'control do not exist.")

    ## test the input files
    seq_info_test = seqinfo(BamFile(test))
    seq_info_control = seqinfo(BamFile(control))
    if(!identical(seq_info_test, seq_info_control))
        warning("The BAM header of 'test' and 'control' differ!")

    suppressWarnings(seqinfo(region) <- seq_info_control[seqlevels(region)])
    region = trim(region)
    chunks = region
    if(max(width(chunks)) > block)
        chunks = subdivideGRanges(chunks, block) ## TODO
    
    header = TRUE
    val = list()
    for(i in seq_along(chunks)) {
        ## the current region
        roi = chunks[i]
        n = width(roi)

        test_tally = tallyBamRegion(test, roi, nCycles, minMapQual)
        control_tally = tallyBamRegion(control, roi, nCycles, minMapQual)
        
        dx = comparativeMismatch(test_tally, control_tally, strand, consensus, roi)

        ## only look at positions with coverage in both samples
        if(select) {
            ## ignore:
            ## 1. no coverage in any sample
            ## 2. NAs
            ## 3. no mismatches in both samples
            idx_good = dx$testDepth > 0 & dx$controlDepth > 0 & !is.na(dx$testDepth) & !is.na(dx$controlDepth) & !is.na(dx$testMismatch) & !is.na(dx$controlMismatch) & (dx$testMismatch > 0 | dx$controlMismatch > 0)
            if(!any(idx_good))
                ## go on if nothing is left
                next
            dx = dx[idx_good, ]
        }
                        
        ## statistics
        pval = with(dx, feTest(testMismatch, testDepth, controlMismatch, controlDepth))
        padj = p.adjust(pval, method = "BH")
        
        ci = with(dx, acCi(testMismatch, testDepth, controlMismatch, controlDepth, beta))

        stat = data.frame(pval = pval, padj = padj)

        ## call significant
        idx_sig = (padj <= alpha) & !is.na(padj)
        idx_ci = ciOutside(ci)
        idx_called = switch(criteria,
            fet = idx_sig,
            ci = idx_ci,
            both = idx_sig & idx_ci,
            any = idx_sig | idx_ci
            )
        ## TODO: How are NAs handeled?
        
        if(select) {
            idx_select = idx_called
        } else {
            idx_select = !logical(n)
        }

        
        dx$called = idx_called
        res = cbind(dx[idx_select, ], ci[idx_select, ], stat[idx_select, ])

        ## classify events
        res = cbind(res, classifyEvent(res))

        ## write output file
        if(write_file && nrow(res) > 0) {
            write.table(res, resultFile, append = !header, quote = FALSE, sep = "\t", row.names = FALSE, col.names = header)
            header = FALSE
        }
        
        if(value) {
            val[[i]] = res
        }        
    }

    res = NULL
    if(value) {
        res = df2gr(do.call(rbind, val))
        seqinfo(res) = seq_info_test[unique(as.character(seqnames(region)))]
        metadata(res) = args
    }
    
    return(res)
}


readRariant <- function(file) {
  
    x = read.table(file, header = TRUE, sep = "\t")
    gr = df2gr(x)
    gr = sort(sortSeqlevels(gr))

    return(gr)
}


writeRariant <- function(x, file) {

    df = gr2df(x)
    
    ## remove and rename columns
    df = df[ ,!(names(df) %in% c("end", "width", "strand"))]
    n = names(df)
    n[match(c("seqnames", "start"), n)] = c("chr", "pos")
    names(df) = n
    
    write.table(df, file, quote = FALSE, sep = "\t", row.names = FALSE)

}


## detect LOH and somatic events ##
classifyEvent <- function(x, alpha = 0.05) {

    pval_null = binomTestPval(x$controlMismatch, x$controlDepth, 0, "greater")
    pval_hetero = binomTestPval(x$controlMismatch, x$controlDepth, 0.5, "both")
    padj_null = p.adjust(pval_null, "BH")
    padj_hetero = p.adjust(pval_hetero, "BH")
    
    is_somatic = (padj_null >= alpha) & (padj_hetero < alpha)
    is_hetero = (padj_hetero >= alpha) & (padj_null < alpha)
    is_powerless = (padj_null >= alpha) & (padj_hetero >= alpha)
          
    verdict = rep("undecided", length(pval_null))
    verdict[is_somatic] = "somatic"
    verdict[is_hetero] = "hetero"
    verdict[is_powerless] = "powerless"
    verdict = factor(verdict, levels = c("somatic", "hetero", "undecided", "powerless"))
    
    res = data.frame(event_type = verdict,
        padj_somatic = padj_null, padj_hetero = padj_hetero,
        pval_somatic = pval_null, pval_hetero = pval_hetero)
    
    return(res)
}


reevalSites <- function(dx, beta = 0.95, alpha = 1 - beta, criteria = c("both", "any", "fet", "ci")) {

    ## statistics
    pval = with(dx, feTest(testMismatch, testDepth, controlMismatch, controlDepth))
    padj = p.adjust(pval, method = "BH")
          
    ci = with(dx, acCi(testMismatch, testDepth, controlMismatch, controlDepth, beta))
  
    stat = data.frame(pval = pval, padj = padj)
    
    ## call significant
    idx_sig = (padj <= alpha) & !is.na(padj)
    idx_ci = ciOutside(ci)
    idx_called = switch(criteria,
        fet = idx_sig,
        ci = idx_ci,
        both = idx_sig & idx_ci,
        any = idx_sig | idx_ci
        )
    ## TODO: How are NAs handeled?
          
    dx$called = idx_called

    ## classify events
    events = classifyEvent(dx)

    ## overwrite existing columns
    new_names = unique(c(names(ci), names(stat), names(events)))
    dxs = dx[ ,!(names(dx) %in% new_names)]
    res = cbind(dxs, ci, stat, events)

    return(res)
}
