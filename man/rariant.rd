\name{rariant}
  
\alias{rariant}
\alias{rariantStandalone}
\alias{readRariant}
\alias{writeRariant}

\title{Rariant calling functions}

\description{
  
  The 'rariant' function screens for variant shifts between a test and
  control sample.  These highlevel functions offers a convenient
  interface for large-scale identification as well as for reexamination
  of existing variant calls.
  
}

\usage{
  
  rariant(test, control, region, beta = 0.95, alpha = 1 - beta,
  select = TRUE, consensus, resultFile, strand = c("both", "plus", "minus"),
  nCycles = 10, minQual = 20, block = 1e4, value = TRUE, criteria = c("both", "any", "fet", "ci"))

  rariantStandalone()
  
  readRariant(file)

  writeRariant(x, file)

}

\arguments{
  
  \item{test, control}{Test and control BAM files.  Other input sources
    will be supported in the future.}
  
  \item{region}{Region(s) of interest to analyze in the calling [GRanges
    with one or more entries].  If missing, the entire genomic space, as
    defined by the BAM headers of the input files, will be covered.}
  
  \item{beta}{Confidence level [numeric in the range [0,1], default:
    0.95].}

  \item{alpha}{Significance threshold for BH-adjusted p-values of the
    Fisher's exact test.}
  
  \item{select}{Should only likely variant positions be selected and
    returned, or the results for all sites be returned.}
  
  \item{consensus}{How to determine the consensus sequence.  By default,
    the consensus is given by the most abundant allele in the control
    sample.  Alternatively, an object with a reference sequence
    ('BSgenome', 'FaFile') can be passed to define the consensus
    sequence.}
  
  \item{resultFile}{If not missing, write the results to a tab-delimited
    file.}
  
  \item{strand}{Which strand should be extracted?  By default, the
    counts of both strands are summed up.}
  
  \item{nCycles}{Number of sequencing cycles to remove from the
    beginning and end of each read when creating the base count
    table. This avoids low quality read positions [default: 10 is
    reasonable for current Illumina sequencing].}
  
  \item{minQual}{Minimum base call quality for reads to be considered
    for the nucleotide count table [default: 20 is reasonable for
    current Illumina sequencing].  Reads with a lower quality are
    dropped.}
  
  \item{block}{Number of the genomic sites to analyze in one chunk. The
    default is a good compromise between memory usage and speed, and
    normally does not require changing.}
  
  \item{value}{Should the results be returned by the function.  For
    calls within R, this is generally set to TRUE and does not need to
    be changed.}

  \item{criteria}{The criteria to determine significant sites.  Criteria
    are: Fisher's exact test, confidence intervals, any or both
    [default] of them.}

  \item{file}{Path to output file from a 'rariant' call.}

  \item{x}{Output of 'rariant' call.}
  
}


\details{

  The 'rariant' function is the workhorse for the comparative variant
  calling and assessment.  It starts with the aligned reads for the test
  (e.g. tumor) and the control (e.g. normal) sample in the BAM format;
  later versions will support additional inputs.

  The 'select' parameter determines whether only significant variant
  sites or all sites are returned.  While the first is suitable for
  detecting variants, the second becomes relevant assessing for example
  the abundance of variants at particular sites of interest - an example
  would be to determine the absence of a specific variant.

  For analyses over large genomic regions and for use with
  infrastructure outside of R, initiating the calling from the command
  line may be a desirable alternative.  The 'rariantStandalone'
  functions returns the full path to a script that can be directly
  called from the command line.  For further details, see the help of
  the script by calling it with the '-h' option, for example 'rariant -h'.

  The 'readRariant' and 'writeRariant' functions allow to import and
  export the results of a 'rariant' call from and to a file output, and
  will return the same object.

}


\value{

  A 'GRanges' object, with each row corresponding to a genomic
  site, and columns:
  
  \itemize{
    
    \item{testMismatch, controlMismatch}{Mismatch counts in the test and
      control sample.}
    
    \item{testDepth, controlDepth}{Sequencing depth in the test and
      control sample.}

    \item{testRef, testAlt}{Reference and alternative allele of the test
      sample.}

    \item{controlRef}{Reference allele of the control sample.}

    \item{testRefDepth, testAltDepth}{Supporting sequencing depth for
      the reference and alternative allele in the test sample.}

    \item{ref}{Consensus allele.}

    \item{p1, p2}{Estimated non-consensus rate in test and control,
      respectively.}

    \item{d}{Estimated shift in the non-consensus rate between test and
      control.}

    \item{ds}{Estimated shift in the non-consensus rate between test and
      control (shrinkage point estimate).}

    \item{lower, upper}{Lower and upper bound of the confidence interval
      for 'd'.}

    \item{pval, padj}{Raw and BH-adjusted p-value of the FET test.}

    \item{called}{Was the site identified as variant?}

    \item{eventType}{The class of the event: somatic, heterozygous,
      undecided.}

    \item{padjSomatic, padjHetero}{BH-adjusted p-values of the
      binomial tests for 'eventType'.}

    \item{pvalSomatic, pvalHetero}{Raw p-values of the binomial tests
      for 'eventType'.}
    
  }
  
}

\examples{

  library(GenomicRanges)
  
  control_bam = system.file("extdata", "NRAS.Control.bam", package = "h5vcData", mustWork = TRUE)
  test_bam = system.file("extdata", "NRAS.AML.bam", package = "h5vcData", mustWork = TRUE)
  
  roi = GRanges("1", IRanges(start = 115258439, end = 115259089))
  
  vars = rariant(test_bam, control_bam, roi)
  
  vars_all = rariant(test_bam, control_bam, roi, select =
  FALSE)

  \dontrun{
    system2(rariantStandalone(), "-h")
  }

}
