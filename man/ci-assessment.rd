\name{ciAssessment}

\alias{coverageProbability}
                
\title{Assessment of CI methods}
 
\description{

  Functions to compute the coverage probability of a confidence interval
  method.

}

\usage{
  coverageProbability(pars, fun = acCi, n_sample = 1e4, min_k, ...)
}

\arguments{

  \item{pars}{Data frame with parameter combinations [data.frame]}
  
  \item{n_sample}{Number of assessments per parameter combination
    [integer(1)].}
  
  \item{fun}{CI function}
  
  \item{min_k}{Minimum 'k2' value to use.}

  \item{...}{Additional arguments that are passed on to 'fun'.}

}

\value{

  The 'data.frame' object 'pars' with additional columns 'cp' for the
  coverage probability and 'aw' average confidence interval width.
  
}

\references{
  
  Fagerland, Morten W., Stian Lydersen, and Petter Laake. Recommended
  Confidence Intervals for Two Independent Binomial Proportions.
  Statistical Methods in Medical Research (2011).
  
}

\examples{
## Define parameter space
pars = expand.grid(k1 = 1:5, k2 = 5, n1 = 30, n2 = 30)
conf_level = 0.95

## Compute coverage probabilities
cp = coverageProbability(pars, fun = acCi, n_sample = 1e2, conf_level = conf_level)
print(cp)
}
