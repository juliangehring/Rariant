\name{ciAdjust}
  
\alias{ciAdjustLevel}
                  
\title{CI Adjust}
   
\description{

  Multiple testing adjustment of confidence levels, as proposed by
  Benjamini and Yekutieli.
  
}
  
\usage{

  ciAdjustLevel(eta0, conf_level)

}
  
\arguments{

  \item{eta0}{Estimated fraction of tests that are consistent with the
    null hypothesis.}

  \item{conf_level}{Unadjusted confidence level}
    
}

\value{

  The adjusted confidence level.

}

\references{

  Benjamini, Yoav, and Daniel Yekutieli. False Discovery Rate–adjusted
  Multiple Confidence Intervals for Selected Parameters. Journal of the
  American Statistical Association 100, no. 469 (2005): 71–81.

}

\examples{
conf_level = 0.95
eta0 = seq(0, 1, by = 0.02)

conf_level_adj = ciAdjustLevel(eta0, conf_level)

plot(eta0, conf_level_adj, pch = 20, ylim = c(conf_level, 1))
}
