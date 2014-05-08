\name{evidenceHeatmap}

\alias{evidenceHeatmap}

\title{Variant Evidence Heatmap}

\description{

  Heatmap with the evidence of variant evidences.
  
}

\usage{

  evidenceHeatmap(x, fill = "d", color = "outside",
                   height = 0.9, width = 0.95, size = 1.5,
                   xvar = "sample", yvar = "loc", ...)

}

\arguments{

  \item{x}{GRanges with variants, as returned by 'rariant'.}

  \item{fill}{Column determining the fill.}

  \item{color}{Column determining the border color.}
  
  \item{height, width}{Height and width of tiles.}
  
  \item{size}{What is this needed for?}
  
  \item{xvar, yvar}{Which column to define the tiles.}
  
  \item{...}{Additional arguments, passed to 'geom_tile'.}
}

\value{

  A ggplot2 object.
  
}

\seealso{

  yesNoMaybe

}
