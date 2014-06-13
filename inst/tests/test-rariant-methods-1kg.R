library(testthat)
library(Rariant)
library(GenomicRanges)
library(ggbio)
library(Rsamtools)

region = GRanges("chr17", IRanges(7565097, 7590856))

## bam files
control_bam = system.file("extdata", "platinum", "control.bam", package = "Rariant", mustWork = TRUE)
test1_bam = system.file("extdata", "platinum", "test.bam", package = "Rariant", mustWork = TRUE)

context("rariant-methods")

v1 = rariant(test1_bam, control_bam, region, select = FALSE)

test_that("'rariant(path, path, GRanges)' works", {
    expect_is(v1, "GRanges")
    expect_identical(length(v1), width(region))
})

context("rariant-methods-2")

v2 = rariant(BamFile(test1_bam), BamFile(control_bam), region, select = FALSE)

context("rariant-methods-2-1")

test_that("'rariant(BamFile, BamFile, GRanges)' works", {
    expect_is(v2, "GRanges")
    expect_identical(length(v2), width(region))
    expect_identical(v1, v2)
})

context("rariant-methods-3")

test_tally = tallyBamRegion(test1_bam, region, 10, 20)
control_tally = tallyBamRegion(control_bam, region, 10, 20)

v3 = rariant(test_tally, control_tally, region, select = FALSE)

context("rariant-methods-3-1")

test_that("'rariant(tally, tally, GRanges)' works", {    
    expect_is(v3, "data.frame")
    expect_identical(nrow(v3), width(region))
})
