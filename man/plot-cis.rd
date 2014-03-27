\name{plotCIs}

\alias{plotCIs}
\alias{plotConfidenceIntervals}
\alias{plotAbundanceShift}
                
\title{Plotting Functions}

\description{

  The 'plotConfidenceIntervals' is a high-level plotting function for
  visualizing confidence intervals.  The 'plotAbundanceShift' function
  visualizes the shift in mismatch rates between two samples.

}

\usage{
plotConfidenceIntervals(x, ylim = c(-1.05, 1.05), color = NULL, ...)

plotAbundanceShift(x, ylim = c(-0.05, 1.05), ...)
}

\arguments{

  \item{x}{'GRanges' with mcols of a CI method, or 'data.frame' as
    returned by one of the CI methods, with the optional column
    'start'.}

  \item{ylim}{Limits of the y-axis.  Using this instead of using the
    'ylim' prevents ugly warnings of 'ggplot2'.}

  \item{color}{Variable that determines the coloring of the confidence
    axis (character).}

  \item{...}{Additional plotting arguments that are passed on to
    ggplot2::geom_pointrange.}

}

\value{

  For a 'GRanges' input: A 'ggbio' object

  For a 'data.frame' input: A 'ggplot' object

}

\examples{
## Generate sample data
counts = data.frame(x1 = 1:5, n1 = 30, x2 = 0:4, n2 = 30)

## Agresti-Caffo
ci_ac = with(counts, acCi(x1, n1, x2, n2))

library(GenomicRanges)
gr = GRanges("1", IRanges(start = 1:nrow(counts), width = 1))
mcols(gr) = ci_ac

## GRanges
plotConfidenceIntervals(gr)

## data.frame
plotConfidenceIntervals(ci_ac)

## abundance shift
plotAbundanceShift(gr)

plotAbundanceShift(ci_ac)
}
