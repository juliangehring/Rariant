comparativeMismatch <- function(test_counts, control_counts, strand = c("both", "plus", "minus"), consensus, roi) {

    #strand = match.arg(strand)

    #test_counts = selectStrand(test_tally, strand)
    #control_counts = selectStrand(control_tally, strand)
    
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
        
    n = nrow(test_counts)
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


dna_bases <- function() {
    res = c("A", "C", "G", "T", "N")
    names(res) = res
    return(res)
}
tallyBamRegion <- function(bam, region, minBase = 0, minMap = 0, maxDepth = 1e4) {

    if(length(region) > 1)
        stop("Region must be of length 1.")

    ## params
    pp = PileupParam(
        max_depth = maxDepth,
        min_base_quality = minBase,
        min_mapq = minMap,
        min_nucleotide_depth = 0,
        min_minor_allele_depth = 0,
        distinguish_strands = FALSE,
        distinguish_nucleotides = TRUE,
        ignore_query_Ns = FALSE,
        include_deletions = FALSE)

    sb = ScanBamParam(which = region)

    t2 = pileup(bam, scanBamParam = sb, pileupParam = pp)

    ## to matrix
    m = matrix(0L, width(region), 4)
    colnames(m) = levels(t2$nucleotide)[1:4]
    idx_col = as.integer(t2$nucleotide)
    if(any(idx_col > 4)) {
        #warning("Skipping non-standard bases")
        t2 = t2[idx_col <= 4, ]
        idx_col = idx_col[idx_col <= 4]
    }
    m[cbind(t2$pos-start(region)+1, idx_col)] = t2$count

    return(m)
}
