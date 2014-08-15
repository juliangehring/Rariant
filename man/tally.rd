\name{tallyBam}

\alias{tallyBamRegion}

\title{Tally a genomic region}

\description{

  Create the nucleotide count table ('tally') of a genomic region from a
  BAM file.

}

\usage{
  tallyBamRegion(bam, region, minBase = 0, minMap = 0, maxDepth = 10000)
}


\arguments{

  \item{bam}{BAM file}

  \item{region}{GRanges with the region to tally, with one entry.}

  \item{minBase, minMap}{Minimum base call and mapping quality for reads
    to be considered for the nucleotide count table [default: 0]. Reads
    with a lower quality are dropped.}

  \item{maxDepth}{Maximal sequencing depth to analyze.}

}


\details{

  For details, look at the documentation of the underlying 'tallyBAM'
  function in the 'h5vc' package.

}


\value{

  An integer array with the dimensions:

  \itemize{

    \item{position}{Length: width(region)}

    \item{base}{A, C, G, T}

    \item{strand}{+, -}

  }
  
}

\seealso{

  h5vc::tallyBAM, deepSNV::bam2R, Rsamtools::pileup
  
}
