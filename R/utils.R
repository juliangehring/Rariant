binGRanges <- function(x, size, method = c("exact", "equal"), trim = FALSE) {

    method = match.arg(method)
    if(method == "exact") {
        ranges = mapply(seq, start(x), end(x), MoreArgs = list(by = size), SIMPLIFY = FALSE)
        names(ranges) = as.character(seqnames(x))
        n = elementLengths(ranges)
        bins = suppressWarnings(
            GRanges(rep(names(ranges), n),
            IRanges(start = unlist(ranges), width = size),
            seqinfo = seqinfo(x))
            )
    } else {
        bins = subdivideGRanges(x, size)
    }
    if(trim)
        bins = trim(bins)

    return(bins)                                                                                                                          
}


grDensity <- function(x, bin_size = 1e+06, method = "exact", chrs = seqlevels(x), seq_gr = as(seqinfo(x), "GRanges")) {

    x = keepSeqlevels(x, chrs)
    bins = binGRanges(seq_gr, bin_size, method)
    bins$counts = countOverlaps(bins, x)
    bins$density = bins$counts/width(bins) * 1e+06

    return(bins)
}
