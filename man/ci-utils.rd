\name{ciUtils}
  
\alias{ciOutside}
\alias{ciCovers}
\alias{ciOverlap}
\alias{ciWidth}
                  
\title{CI Utils}
   
\description{

  Utility functions to find confidence intervals that (a) overlap a
  certain value ('ciOutside', 'ciCovers') and (b) different confidence
  intervals overlap ('ciOverlap').
  
}
  
\usage{

  ciOutside(x, delta = 0)

  ciCovers(x, delta = 0)

  ciOverlap(x, y)

  ciWidth(x)
}
  
\arguments{
  
  \item{x, y}{CIs, as obtained from e.g. the 'acCi' function.}  

  \item{delta}{Variant frequency value to check against [default: 0].}
    
}

\value{

  A logical vector, where each elements corresponds to the respective
  row of 'x' (and 'y').

  For 'ciWidth': A numeric vector with the widths of the confidence
  intervals.

}

\examples{
## Generate sample data
counts = data.frame(x1 = 1:5, n1 = 30, x2 = 0:4, n2 = 30)
  
## Agresti-Caffo
ci_ac = with(counts, acCi(x1, n1, x2, n2))
ci_ac2 = with(counts, acCi(x1, n1, x2, n2, 0.99))
  
## cover 0
idx_zero = ciCovers(ci_ac)

## cover 1
idx_one = ciCovers(ci_ac, delta = 1)

## overlap
idx_same = ciOverlap(ci_ac, ci_ac2)

## width
width = ciWidth(ci_ac)
}
