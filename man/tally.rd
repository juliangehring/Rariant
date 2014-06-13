\name{tallyBam}

\alias{tallyBamRegion}
\alias{tallyBamPart}              

\title{Tally a genomic region}

\description{

  Create the nucleotide count table ('tally') of a genomic region from a
  BAM file.

}

\usage{
  tallyBamRegion(file, region, nCycles = 0, minQual = 0)

  tallyBamPart(file, chrom, start, end, nCycles = 0, minQual = 0)
}


\arguments{

  \item{file}{BAM file path}

  \item{region}{GRanges with the region to tally, with one entry.}

  \item{chrom, start, end}{Alternative to a 'GRanges' input, the
    chromosome, start, and end positon.}

  \item{nCycles}{Number of sequencing cycles to remove from the
    beginning and end of each read when creating the base count
    table. This avoids low quality read positions [default: 0].}
  
  \item{minQual}{Minimum base call quality for reads to be considered
    for the nucleotide count table [default: 0]. Reads with a lower
    quality are dropped.}

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

  h5vc::tallyBAM, deepSNV::bam2R
  
}
