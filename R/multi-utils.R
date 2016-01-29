## what is it good for?
.combineCalls <- function(...) {
    l = GenomicRangesList(...)
    sample_names = names(l)
    if(is.null(sample_names))
        sample_names = as.character(seq_along(l))
    return(l)
}

mergeCalls <- function(x) {
    sample_names = names(x)
    if(is.null(sample_names))
        sample_names = as.character(seq_along(x))
    f <- function(x) {
        x1 = as(x, "data.frame")
        x1$n = 1:nrow(x1)
        x1$outside = ciOutside(x1)
        return(x1)
    }
    lx = lapply(x, f)
    z = do.call(rbind, lx)
    #z = rbind_all(lx)
    ## Has problems if corresponding colums have different types
    z$sample = factor(rep(sample_names, elementNROWS(x)),
        levels = sample_names)
    return(z)
}

updateCalls <- function(x, ...) {
    x$verdict = yesNoMaybe(x, ...)
    x$loc = factor(seq_along(x))
    return(x)
}

filterCalls <- function(x, ..., minCount = 1) {
    ind = findCalls(x, ..., minCount = minCount)
    y = endoapply(x, "[", ind)
    return(y)
}

findCalls <- function(x, ..., minCount = 1) {
    idx_tab= sapply(x, function(y, ...) {y %in% subset(y, ...)}, ...)
    ind = which(rowSums(idx_tab) >= minCount)
    return(ind)
}

checkCalls <- function(x) {
    ns = elementNROWS(x)
    stopifnot(length(unique(ns)) == 1)
}
