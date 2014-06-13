exampleCalls <- function(select = FALSE) {

    tp53_region = GRanges("chr17", IRanges(7565097, 7590856))

    control_bam = system.file("extdata", "platinum", "control.bam", package = "Rariant", mustWork = TRUE)
    test1_bam = system.file("extdata", "platinum", "test.bam", package = "Rariant", mustWork = TRUE)
    test2_bam = system.file("extdata", "platinum", "test2.bam", package = "Rariant", mustWork = TRUE)
    mix_bam = system.file("extdata", "platinum", "mix.bam", package = "Rariant", mustWork = TRUE)

    v_test1 = rariant(test1_bam, control_bam, tp53_region, select = select)
    v_test2 = rariant(test2_bam, control_bam, tp53_region, select = select)
    v_mix = rariant(mix_bam, control_bam, tp53_region, select = select)

    l = GenomicRangesList(t1 = v_test1, t2 = v_test2, m = v_mix)

    return(l)
}
