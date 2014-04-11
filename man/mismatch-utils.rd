\name{mismatchUtils}
      
\alias{mismatchUtils}
\alias{selectStrand}
\alias{seqDepth}
\alias{callConsensus}
\alias{mismatchCount}

\title{Tally processing low-level functions}

\description{
  
  Functions for processing position-specific base count tables (tallies)
  and extracting mismatches counts.
  
}
      
\usage{

## low-level functions
selectStrand(x, strand = c("both", "plus", "minus"), idx = 1:ncol(x))
    
seqDepth(x)
    
callConsensus(counts, verbose = FALSE)
    
mismatchCount(counts, consensus, depth = rowSums(counts))

}
      
\arguments{

  \item{x}{Input object}

  \item{strand}{Which strand to return?}

  \item{idx}{Index of bases to consider (leave as is)}

  \item{counts}{Count matrix}

  \item{verbose}{Show warnings}

  \item{consensus}{Consensus sequence}

  \item{depth}{Sequencing depth for counts}.
  
}

\seealso{

  comparativeMismatch
  
}
