\name{splitSample}

\alias{splitSample}
\alias{splitSampleBinom}
                
\title{Split Sample for Binomial Data}

\description{Sample splitting, according to Hall, 2014.}

\usage{
splitSampleBinom(x, n)
}

\arguments{

  \item{x}{Number of successes}

  \item{n}{Number of trials}

}

\details{
  
  These functions implement sample splitting of a binomial rate.

  Note that the results depend on the state of the random number
  generator, and are therefore not strictly deterministic.
  
}

\value{
  
  A vector with the rate \eqn{p = \frac{X}{N}}{p = X/N}, obtained with
  sample splitting.
  
}

\references{

  Decrouez, Geoffrey, and Peter Hall. "Split Sample Methods for
  Constructing Confidence Intervals for Binomial and Poisson
  Parameters." Journal of the Royal Statistical Society: Series B
  (Statistical Methodology), 2013, n/a–n/a. doi:10.1111/rssb.12051.

}

\examples{
n = 10
m = 5
pt = 0.5

x = rbinom(m, n, pt)
p = x/n

ps = splitSampleBinom(x, n)

round(cbind(p, ps), 2)
}
