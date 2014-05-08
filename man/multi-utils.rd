\name{multiCalls}

\alias{findCalls}
\alias{filterCalls}
\alias{mergeCalls}
\alias{updateCalls}

\title{Multi call processing}

\description{

  Utilities for processing matched calls from multipe samples.
  
}

\usage{

  findCalls(x, ..., minCount = 1)
  filterCalls(x, ..., minCount = 1)

  mergeCalls(x)

  updateCalls(x, ...)
}


\arguments{

  \item{x}{GenomicRangesList with calls from multiple samples.}
  
  \item{...}{Additional arguments.}

  \item{minCount}{For finding and filtering, for how many samples must
    the condition '...' hold true for a site to be returned?}

}
