\name{yesNoMaybe}

\alias{yesNoMaybe}

\title{Determine Variant Evidence}

\description{

  Determine the evidence (absence, presence, dontknow) of variants.

}

\usage{
  yesNoMaybe(x, null = 0, one = 0.5)
}


\arguments{

  \item{x}{GRanges with variants, as returned by 'rariant'.}

  \item{null}{Shift consistent with the _absence_ of a variant.}

  \item{one}{Shift consistent with the _presence_ of a variant.}
  
}

\value{

  The same GRanges object as the input 'x', with the factor column
  'verdict': 'absent', 'present', 'inbetween', 'dontknow'
  
}
