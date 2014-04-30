library(testthat)
library(Rariant)
library(GenomicRanges)

context("rariant")

control_bam = system.file("extdata", "NRAS.Control.bam", package = "h5vcData", mustWork = TRUE)
test_bam = system.file("extdata", "NRAS.AML.bam", package = "h5vcData", mustWork = TRUE)
    
roi = GRanges("1", IRanges(start = 115258439, end = 115259089))

vars = rariant(test_bam, control_bam, roi)

vars_all = rariant(test_bam, control_bam, roi, select = FALSE)

test_that("'rariant' works", {

    expect_is(vars, "GRanges")
    expect_is(vars_all, "GRanges")

    expect_equal(length(vars), 1L)
    expect_equal(length(vars_all), width(roi))

    expect_error(rariant(test = test_bam, region = roi))
    expect_error(rariant(control = control_bam, region = roi))
    expect_error(rariant())
    
})

context("plots")  

test_that("'plotConfidenceIntervals' works", {

    p = plotConfidenceIntervals(vars)
    p_col = plotConfidenceIntervals(vars, color = "eventType")

    expect_is(p, "GGbio")
    expect_is(p_col, "GGbio")

})


context("IO")

test_that("'readRariant/writeRariant' works", {

    tf = tempfile()
    writeRariant(vars_all, tf)
    r = readRariant(tf)

    expect_equal(dim(mcols(vars_all)), dim(mcols(r)))
    expect_equal(names(mcols(vars_all)), names(mcols(r)))
    
})
