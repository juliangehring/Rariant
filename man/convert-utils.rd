\name{convertUtils}
  
\alias{gr2pos}
\alias{pos2gr}
                  
\title{Position converters}
   
\description{

  Utility functions to convert between 'GRanges' and 'character' objects. 
  
}
  
\usage{

gr2pos(x)

pos2gr(x)

}
  
\arguments{
  
  \item{x}{GRanges or character object.}  
    
}

\value{

  A GRanges object or character object, with the position.

}

\examples{
library(GenomicRanges)

gr = GRanges(1:2, IRanges(1:2, width = 1))

pos = gr2pos(gr)
gr2 = pos2gr(pos)

identical(gr, gr2)
}
