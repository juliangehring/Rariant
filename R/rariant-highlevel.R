readRariant <- function(file, ...) {
  
    x = read.table(file, header = TRUE, sep = "\t", ...)
    gr = df2gr(x)
    gr = sort(sortSeqlevels(gr))

    return(gr)
}


writeRariant <- function(x, file, ...) {

    df = gr2df(x)
    
    ## remove and rename columns
    df = df[ ,!(names(df) %in% c("end", "width", "strand"))]
    n = names(df)
    n[match(c("seqnames", "start"), n)] = c("chr", "pos")
    names(df) = n
    
    write.table(df, file, quote = FALSE, sep = "\t", row.names = FALSE, ...)

}


## detect LOH and somatic events ##
classifyEvent <- function(x, alpha = 1e-3, offset = 0.5) {

    pval_null = binomTestPval(pmin(x$controlMismatch, x$controlDepth-x$controlMismatch),
        x$controlDepth, offset/x$controlDepth, "greater")
    pvalHetero = binomTestPval(x$controlMismatch, x$controlDepth, 0.5, "both")
    padj_null = p.adjust(pval_null, "BH")
    padjHetero = p.adjust(pvalHetero, "BH")
    
    is_somatic = (pval_null >= alpha) & (pvalHetero < alpha)
    is_hetero = (pvalHetero >= alpha) & (pval_null < alpha)
    is_powerless = (pval_null >= alpha) & (pvalHetero >= alpha)
          
    verdict = rep("undecided", length(pval_null))
    verdict[is_somatic] = "somatic"
    verdict[is_hetero] = "hetero"
    verdict[is_powerless] = "powerless"
    verdict = factor(verdict, levels = c("somatic", "hetero", "undecided", "powerless"))
    
    res = data.frame(eventType = verdict,
        padjSomatic = padj_null, padjHetero = padjHetero,
        pvalSomatic = pval_null, pvalHetero = pvalHetero)
    
    return(res)
}


reevalSites <- function(dx, beta = 0.95, alpha = 1 - beta, criteria = c("both", "any", "fet", "ci")) {

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
