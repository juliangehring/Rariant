\name{mismatchPlotRegion}

\alias{mismatchPlotRegion}
\alias{tallyPlot}

\title{Mismatch plot from BAM files}

\description{

  Create a mismatch plot from a list of BAM files directly.

}

\usage{
  mismatchPlotRegion(files, region, win = 25, ref, ...)

  tallyPlot(file, region, ref, nCycles = 0, minQual = 0, minFreq = 0, ...)
}


\arguments{

  \item{files, file}{BAM file paths}

  \item{region}{GRanges with the position (width: 1) to tally, with one
    entry.}

  \item{win}{Window size, extending the region on both sides.}
  
  \item{ref}{Reference object, as 'BSgenome'.}

  \item{nCycles}{Number of sequencing cycles to remove from the
    beginning and end of each read when creating the base count
    table. This avoids low quality read positions [default: 0].  See
    'tallyBamRegion'}
  
  \item{minQual}{Minimum base call quality for reads to be considered
    for the nucleotide count table [default: 0]. Reads with a lower
    quality are dropped.  See 'tallyBamRegion'}

  \item{minFreq}{Currently not used}
  
  \item{...}{Additional arguments, passed to 'tallyBAM'.}
  
}


\value{

  A 'ggplot2' or 'ggbio' object.
  
}

\seealso{

  h5vc::mismatchPlot
  
}

\examples{
  library(ggbio)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)

  region = GRanges("chr17", IRanges(7572100, width = 1))

  control_bam = system.file("extdata", "control.bam", package =
  "Rariant", mustWork = TRUE)
  mix_bam = system.file("extdata", "mix.bam", package = "Rariant",
  mustWork = TRUE)

  bam_files = c(control_bam, mix_bam)

  p = mismatchPlotRegion(bam_files, region, win = 50,
  BSgenome.Hsapiens.UCSC.hg19)

  print(p)

  region = GRanges("chr17", IRanges(7572050, width = 100))

  control_bam = system.file("extdata", "control.bam", package =
    "Rariant", mustWork = TRUE)
  test1_bam = system.file("extdata", "test.bam", package =
    "Rariant", mustWork = TRUE)
  test2_bam = system.file("extdata", "test2.bam", package =
    "Rariant", mustWork = TRUE)
  mix_bam = system.file("extdata", "mix.bam", package =
    "Rariant", mustWork = TRUE)

  bam_files = c(control_bam, test1_bam, test2_bam, mix_bam)

  library(BSgenome.Hsapiens.UCSC.hg19)
  ref = BSgenome.Hsapiens.UCSC.hg19

  p = tracks(lapply(bam_files, tallyPlot, region, ref, minQual = 25))

  print(p)
}
