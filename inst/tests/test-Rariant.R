library(testthat)
library(Rariant)

context("acCi")

test_that("'acCi' works", {

    ## from the Agresti paper
    test_df = data.frame(
        x1 = c(5, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2),
        x2 = c(0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4),
        n = 10,
        lower = c(0.093, -0.020, -0.124, -0.222, -0.314, -0.400, -0.124, -0.240, -0.346, -0.446, -0.538),
        upper = c(0.740, 0.686, 0.624, 0.556, 0.481, 0.400, 0.457, 0.407, 0.346, 0.279, 0.205))

    ## CIs
    ci_est = with(test_df, acCi(x1, n, x2, n, 0.95))

    expect_equal(ci_est[ ,"lower"], test_df[ ,"lower"], tolerance = 0.01)
    expect_equal(ci_est[ ,"upper"], test_df[ ,"upper"], tolerance = 0.01)

})


test_that("'acCi' returns 'NA' for sample sizes of 0", {
 
    test_df = data.frame(
        x1 = c(0, 1, 1, 1),
        x2 = c(0, 1, 1, 1),
        n1 = c(0, 0, 0, 1),
        n2 = c(0, 0, 1, 0)
        )
    
    ## CIs
    ci_est = with(test_df, acCi(x1, n1, x2, n2, 0.95))
    
    ixd = all(is.na(ci_est[ ,"lower"])) & all(is.na(ci_est[ ,"upper"]))
    
    expect_true(ixd)
    
})


test_that("'nhsCi' returns 'NA' for sample sizes of 0", {
     
    test_df = data.frame(
        x1 = c(0, 1, 1, 1),
        x2 = c(0, 1, 1, 1),
        n1 = c(0, 0, 0, 1),
        n2 = c(0, 0, 1, 0)
        )

    ## CIs
    ci_est = with(test_df, nhsCi(x1, n1, x2, n2, 0.95))
    
    ixd = all(is.na(ci_est[ ,"lower"])) & all(is.na(ci_est[ ,"upper"]))
    
    expect_true(ixd)
    
})


context("rariant")

library(GenomicRanges)

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
