\name{rariantInspect}
  
\alias{rariantInspect}

\title{Interactive inspection}

\description{
  
  Interactively inspect variant sites and results of the 'rariant'
  function.

}

\usage{
    
  rariantInspect(x)

}

\arguments{

  \item{x}{The return value of the 'rariant' or 'readRariant' function.}
  
}
  

\details{

  With the web interface of 'rariantInspect' can existing variant calls
  and assessment be explored interactively.  It allows to select the
  genomic region of interest and the type of event.  Results are shown
  as both a confidence interval plot and a results table that can be
  further filtered and reordered.

}


\examples{
\donttest{
example(rariant)

rariantInspect(vars_all)
}
}
