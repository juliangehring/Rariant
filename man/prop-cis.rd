\name{propCIs}

\alias{propCIs}
\alias{acCi}
\alias{nhsCi}
                
\title{Confidence Interval Functions}

\description{Vectorized implementation of confidence intervals}

\usage{
acCi(x1, n1, x2, n2, conf_level = 0.95, clip = TRUE, split = FALSE)
nhsCi(x1, n1, x2, n2, conf_level = 0.95)
}

\arguments{

  \item{x1}{Mismatch counts in the test sample.}

  \item{n1}{Sequencing depth (total counts) in the test sample.}

  \item{x2}{Mismatch counts in the control sample.}

  \item{n2}{Sequencing depth (total counts) in the control sample.}

  \item{conf_level}{Confidence level $beta$ (default: 0.95).}

  \item{clip}{Should the CIs be clipped to the interval [-1,1] if they
    exceed this?}

  \item{split}{Should the sample split method be applied?  See
    'splitSampleBinom' for details.}
  
}

\details{
  
  These functions implement a vectorized version of the two-sided
  Agresti-Caffo, and Newcombe-Hybrid-Score confidence interval for the
  difference of two binomial proportions.

}

\value{

  A data frame with columns

  \itemize{

    \item{d}{Estimate for the difference of rates 'p1' and 'p2'.}

    \item{p1, p2}{Estimates for the mismatches rates for each sample.}

    \item{lower, upper}{Lower and upper bound of the confidence
      interval.}

    \item{w}{Width of the confidence interval.}

  }
  
}

\references{

  Agresti, Alan, and Brian Caffo. Simple and Effective Confidence
  Intervals for Proportions and Differences of Proportions Result from
  Adding Two Successes and Two Failures. The American Statistician 54,
  no. 4 (2000): 280–288

  Newcombe, Robert G. Interval Estimation for the Difference between
  Independent Proportions: Comparison of Eleven Methods. Statistics in
  Medicine 17, no. 8 (1998): 873–890.

  Fagerland, Morten W., Stian Lydersen, and Petter Laake. Recommended
  Confidence Intervals for Two Independent Binomial Proportions.
  Statistical Methods in Medical Research (2011).

}

\seealso{

  nhsCi

  splitSampleBinom
  
  binMto::Add4
  binMto::NHS
  
}

\examples{
## Generate sample data
counts = data.frame(x1 = 1:5, n1 = 30, x2 = 0:4, n2 = 30)

## Agresti-Caffo
ci_ac = with(counts, acCi(x1, n1, x2, n2))

## Newcombe-Hybrid Score
ci_nhs = with(counts, nhsCi(x1, n1, x2, n2))

print(ci_ac)
}
