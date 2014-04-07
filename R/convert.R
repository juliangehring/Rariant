df2gr <- function(x) {
    n = c("chr", "pos")
    vr = GRanges(x$chr, IRanges(x$pos, width = 1L))
    mcols(vr) = x[ ,!names(x) %in% n]
    #    ref = x$test_ref, alt = x$test_alt,
    #    refDepth = x$ref_depth, altDepth = x$alt_depth, totalDepth = x$test_depth)
    #n = c("chr", "pos", "test_ref", "test_alt", "ref_depth", "alt_depth")
    #mcols(vr) = x[ ,!(names(x) %in% n)]
    return(vr)
}


gr2df <- function(x) {
    res = as(x, "data.frame")
    return(res)
}


gr2vr <- function(gr) {
    vr = VRanges(seqnames(gr), ranges(gr))
    mcols(vr) = mcols(gr)
    return(vr)
}


vr2gr <- function(vr) {
    gr = as(vr, "GRanges")
    return(gr)
}


gr2pos <- function(x) {
    res = sprintf("%s:%d-%d", as.character(seqnames(x)), start(x), end(x))
    return(res)
}


pos2gr <- function(x) {
    ## this is optimized
    y = unlist(strsplit(x, ";"))
    idx_range = grepl("-", y) ## ranges have a '-'
    ## match ranges
    if(any(idx_range)) {
        m = str_capture_matrix(y[idx_range], "(.+):(\\d+)-(\\d+)", 3)
        res1 = GRanges(m[ ,2], IRanges(as.integer(m[ ,3]), as.integer(m[ ,4])))
    } else {
        res1 = GRanges()
    }
    ## match positions
    if(any(!idx_range)) {
        m = str_capture_matrix(y[!idx_range], "(.+):(\\d+)", 2)
        res2 = GRanges(m[ ,2], IRanges(as.integer(m[ ,3]), as.integer(m[ ,3])))
    } else {
        res2 = GRanges()
    }
    ## restore original order, fast enough
    ord = c(which(idx_range), which(!idx_range))
    res = c(res1, res2)[ord]
    return(res)
}

## faster than str_match
str_capture_matrix <- function(x, reg, n) {
    r = regexec(reg, x)
    l = regmatches(x, r)
    g = matrix(unlist(l), ncol = n + 1, byrow = TRUE)
    return(g)
}

  
ds2int <- function(ds) {
    ind = as.numeric(ds)
    ind[ind == 15] = 16
    ind = as.integer(log2(ind) + 1)
    return(ind)
}

  
mat2ind <- function(i, n) {
    return( (i - 1) * n + 1:n )
}


ind2base <- function(ind, bases = dna_bases()) {
    ind[is.na(ind)] = 5
    res = factor(bases[ind], levels = bases)
    return(res)
}
