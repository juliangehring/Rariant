rariantFromMatrix <- function(
    test,
    control,
    region,
    beta = 0.95,
    alpha = 1 - beta,
    select = TRUE,
    consensus,
    strand = c("both", "plus", "minus"),
    criteria = c("both", "any", "fet", "ci")
    ) {
    
    ## input arguments
    criteria = match.arg(criteria)
    roi = region

    if(!identical(dim(test), dim(control)))
        warning("The BAM headers of 'test' and 'control' differ!")

    n = dim(test)[1]

    test_tally = test
    control_tally = control
    
    dx = comparativeMismatch(test_tally, control_tally, strand, consensus, roi)
    
    ## only look at positions with coverage in both samples
    if(select) {
        ## ignore:
        ## 1. no coverage in any sample
        ## 2. NAs
        ## 3. No mismatches in both samples (no chance for a change)
        ## not optimal in terms of speed, but okay
        idx_depth = dx$testDepth > 0 & dx$controlDepth > 0 & !is.na(dx$testDepth) & !is.na(dx$controlDepth)
        idx_mismatch = !is.na(dx$testMismatch) & !is.na(dx$controlMismatch) & (dx$testMismatch > 0 | dx$contolMismatch > 0)
        idx_good = idx_depth & idx_mismatch
        if(!any(idx_good))
            ## go on if nothing is left
            return(NULL)
        dx = dx[idx_good, ]
    }
    
    ## statistics
    pval = feTest(dx$testMismatch, dx$testDepth, dx$controlMismatch, dx$controlDepth)
    padj = p.adjust(pval, method = "BH")
    
    ci = acCi(dx$testMismatch, dx$testDepth, dx$controlMismatch, dx$controlDepth, beta)
    
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
    
    if(select) {
        idx_select = idx_called
    } else {
        idx_select = !logical(n)
    }
    
    dx$called = idx_called
    res = cbind(dx[idx_select, ], ci[idx_select, ], stat[idx_select, ])
    
    ## classify events
    res = cbind(res, classifyEvent(res))

    return(res)
}
