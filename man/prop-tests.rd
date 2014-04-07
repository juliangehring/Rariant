\name{propTests}

\alias{propTests}
\alias{scoreTest}
\alias{nmTest}
\alias{feTest}
                
\title{Testing Functions}

\description{

  Vectorized implementation of testing functions

}

\usage{

scoreTest(x1, n1, x2, n2)
  
nmTest(x1, n1, x2, n2, delta = 0)

feTest(x1, n1, x2, n2, ...)
}

\arguments{

  \item{x1}{Mismatch counts in the test sample.}

  \item{n1}{Sequencing depth (total counts) in the test sample.}

  \item{x2}{Mismatch counts in the control sample.}

  \item{n2}{Sequencing depth (total counts) in the control sample.}

  \item{delta}{Difference to test against (default: 0).}

  \item{...}{Additional arguments.}
}

\details{
  
  These functions implement a vectorized version of the two-sided (a)
  Score test and (b) Miettinen-Nurminen test for the difference between
  to Binomial proportions.

  Usage of the score test is discouraged in the settings considered
  here, since it is ill-defined for positions with no mismatches.

}

\value{

  A data frame with columns

  \itemize{

    \item{dhat}{Estimate for the difference of rates 'p1' and 'p2'.}

    \item{p1, p2}{Estimates for the mismatches rates for each sample.}

    \item{tval}{T-value}

    \item{pval}{P-value}

  }
  
}

\references{

  Miettinen, Olli, and Markku Nurminen. Comparative Analysis of Two
  Rates. Statistics in Medicine 4, no. 2 (1985):
  213–226. doi:10.1002/sim.4780040211.

}

\seealso{

  VariantTools package
  
}

\examples{
## Generate sample data
counts = data.frame(x1 = 1:5, n1 = 30, x2 = 0:4, n2 = 30)

## Score test
stat_st = with(counts, scoreTest(x1, n1, x2, n2))

## NM test
stat_nm = with(counts, nmTest(x1, n1, x2, n2))

## Fisher test
stat_fet = with(counts, feTest(x1, n1, x2, n2))

print(stat_st)

print(stat_nm)

print(stat_fet)
}
