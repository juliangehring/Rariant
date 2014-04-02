comparativeMismatch <- function(test_tally, control_tally, strand = c("both", "plus", "minus"), consensus, roi) {

    strand = match.arg(strand)

    test_counts = selectStrand(test_tally, strand)
    control_counts = selectStrand(control_tally, strand)
    
    testDepth = seqDepth(test_counts)
    controlDepth = seqDepth(control_counts)
    
    test_base = callConsensus(test_counts) ## needed?
    control_base = callConsensus(control_counts)

    ## define the consensus sequence
    if(!missing(consensus)) {
        ref_base = suppressWarnings(ds2int(getSeq(consensus, roi)[[1]]))
        ref_base[ref_base == 5] = NA ## 'N' -> NA
    } else {
        ref_base = control_base
    }
     
    ## mismatch counts
    controlMismatch = mismatchCount(control_counts, ref_base, controlDepth)
    controlMismatch[is.na(controlMismatch)] = 0 ## TODO
    testMismatch = mismatchCount(test_counts, ref_base, testDepth)
    testMismatch[is.na(testMismatch)] = 0
        
    n = nrow(test_tally)
    ## depth of the consensus allele in test
    refDepth = test_counts[mat2ind(ref_base, n)]
    refDepth[is.na(refDepth)] = 0

    ## substract the consensus counts for test
    idx_consensus = mat2ind(ref_base, n)
    test_counts[idx_consensus] = 0 ## [cbind(i, j)] as alternative
        
    ## call alt alleles and counts
    test_alt_base = callConsensus(test_counts)
    altDepth = test_counts[mat2ind(test_alt_base, n)]

    ## TODO: make this more efficient
    dx = data.frame(
        chr = seqchar(roi), pos = start(roi):end(roi),
        testMismatch = testMismatch, controlMismatch = controlMismatch,
        testDepth = testDepth, controlDepth = controlDepth,
        testRef = ind2base(test_base), testAlt = ind2base(test_alt_base),
        controlRef = ind2base(control_base), ref = ind2base(ref_base),
        refDepth = refDepth, altDepth = altDepth
        )

    return(dx)
}

  
selectStrand <- function(x, strand = c("both", "plus", "minus"), idx = 1:ncol(x)) {
    if(strand == "both")
        y = x[ ,idx,1] + x[ ,idx,2]
    if(strand == "plus")
        y = x[ ,idx,1]
    if(strand == "minus")
        y = x[ ,idx,2]
    dim(y) = c(nrow(x), length(idx))
    return(y)
}


seqDepth <- function(x) {
      return(as.integer(rowSums(x)))
}


callConsensus <- function(counts, verbose = FALSE) {
    idx_max = max.col(counts, ties.method = "first")
    idx_max_2 = max.col(counts, ties.method = "last")
    idx = (idx_max == idx_max_2)
    if(verbose && !all(idx))
        warning(sprintf("%d %s", sum(!idx), "calls are ambiguous."))
    idx_max[!idx] = NA
    return(idx_max)
}


mismatchCount <- function(counts, consensus, depth = rowSums(counts)) {
    idx_mat = (consensus - 1) * nrow(counts) + 1:nrow(counts)
    mmc = (depth - counts[idx_mat])
    return(mmc)
}


#resolveConsensus <- function(..., method = c("ref", "min")) {    
#    method = match.arg(method)
#    idx_amb = is.na(control_base)
#    if(any(idx_amb)) {
#        if(method == "min") {        
#            ter = sapply(seq_along(idx_bases),
#                function(i, r, c) {mismatchCount(r, i, c)},
#                test_counts[idx_amb, ], test_cov[idx_amb])
#            control_base[idx_amb] = callConsensus(-ter)
#        }
#        if(method == "ref") {
#            control_base[idx_amb] = ref_base[idx_amb]
#        }
#    }
#}

dna_bases <- function() {
    res = c("A", "C", "G", "T", "N")
    names(res) = res
    return(res)
}


tallyBamRegion <- function(file, region, ncycles, minq) {

    tally = tallyBAM(file,
            chr = as.character(seqnames(region)),
            start = start(region), stop = end(region),
            q = minq, ncycles = ncycles)
    ## only HQ counts
    tally = aperm(tally[5:8, , ,drop=FALSE], c(3, 1, 2)) ## dont drop if only one pos

    return(tally)
}
