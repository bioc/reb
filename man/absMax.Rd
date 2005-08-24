\name{absMax}
\alias{absMax}
\title{ Absolute Maxima}
\description{
  Returns the absolute maxima of the input values. 
}
\usage{
absMax(x)
}
\arguments{
  \item{x}{numeric argument}
}
\value{
  \code{absMax} returns the absolute maximum of all the values present
  in the arguments as \code{double} persurving the sign. Essentially
  \code{max(abs(x),na.rm=T)}.
}
\author{Karl A. Dykema and Kyle A. Furge}
\examples{
  absMax(c(1,2,3,4))
  absMax(c(-1,-2,-3,-4))
}
\keyword{ manip }
