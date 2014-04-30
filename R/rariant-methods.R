setGeneric("rariant",
           function(test, control, region, ...)
           standardGeneric("rariant")
           )

setMethod("rariant",
          signature = signature(test = "BamFile", control = "BamFile", region = "GRanges"),
          function(test, control, region, beta = 0.95, alpha = 1 - beta, select = TRUE, consensus, resultFile, strand = c("both", "plus", "minus"), nCycles = 10, minQual = 20, block = 1e4, value = TRUE, criteria = c("both", "any", "fet", "ci")) {
              rariantFromBam(test, control, region, beta, alpha, select, consensus, resultFile, strand, nCycles, minQual, block, value, criteria)
          })

setMethod("rariant",
          signature = signature(test = "character", control = "character", region = "GRanges"),
          function(test, control, region, beta = 0.95, alpha = 1 - beta, select = TRUE, consensus, resultFile, strand = c("both", "plus", "minus"), nCycles = 10, minQual = 20, block = 1e4, value = TRUE, criteria = c("both", "any", "fet", "ci")) {
              rariantFromBam(BamFile(test), BamFile(control), region, beta, alpha, select, consensus, resultFile, strand, nCycles, minQual, block, value, criteria)
          })

setMethod("rariant",
          signature = signature(test = "array", control = "array", region = "GRanges"),
          function(test, control, region, beta = 0.95, alpha = 1 - beta, select = TRUE, consensus, strand = c("both", "plus", "minus"), criteria = c("both", "any", "fet", "ci")) {
              rariantFromMatrix(test, control, region, beta, alpha, select, consensus, strand, criteria)
          })
