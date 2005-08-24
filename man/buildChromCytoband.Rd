\name{buildChromCytoband}
\alias{buildChromCytoband}

\title{ Construct a chromLocation object from a cytoband environment}
\description{
Construct a chromLocation object from a cytoband environment. Human, Rat, and Mouse are currently possible.
}
\usage{
buildChromCytoband(organism = "h")
}

\arguments{
  \item{organism}{ character, "h" for human, "m" for mouse, and "r" for rat. }
}

\value{
a chromLocation object
}
\author{ 
  Karl J. Dykema, \email{karl.dykema@vai.org} 
  Kyle A. Furge, \email{kyle.furge@vai.org}
}

\seealso{ ~\code{\link{buildChromLocation}} }
\examples{
  humanBands <- buildChromCytoband("h")
  humanBands@chromLocs[["1"]]
  
}
\keyword{ manip}
