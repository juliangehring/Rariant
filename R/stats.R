call_base <- function(counts, verbose = FALSE) {
    idx_max = max.col(counts, ties.method = "first")
    idx_max_2 = max.col(counts, ties.method = "last")
    idx = (idx_max == idx_max_2)
    if(verbose && !all(idx))
        warning(sprintf("%d %s", sum(!idx), "calls are ambiguous."))
    idx_max[!idx] = NA
    return(idx_max)
}

binomTestPval <- function(x, n, p, alternative = c("less", "greater", "both")) {

    alternative = match.arg(alternative)
    if(alternative == "both" && p != 0.5)
        warning("Approximated two sided Binomial test")
    pval = switch(alternative,
        less = pbinom(x, n, p),
        greater = pbinom(x - 1, n, p, lower.tail = FALSE),
        both = pmin(2*pmin(pbinom(x, n, p), pbinom(x - 1, n, p, lower.tail = FALSE)), 1),
        )

    return(pval)
}
