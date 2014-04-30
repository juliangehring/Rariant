library(testthat)
library(Rariant)
library(GenomicRanges)

context("tally")

test_that("'tallyBamRegion' works", {
    
    #control_bam = system.file("extdata", "NRAS.Control.bam", package = "h5vcData", mustWork = TRUE)
    test_bam = system.file("extdata", "NRAS.AML.bam", package = "h5vcData", mustWork = TRUE)
    
    roi = GRanges("1", IRanges(start = 115258439, end = 115259089))

    test_tally = tallyBamRegion(test_bam, roi)
    expect_is(test_tally, "array")
    expect_identical(dim(test_tally), c(width(roi), 4L, 2L))
    expect_identical(names(dimnames(test_tally)), c("position", "base", "strand"))

})
