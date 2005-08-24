\name{naMean}
\alias{naMean}
\title{Wrapper function for the arithmetic mean}
\description{
  Simple call to mean with the \code{na.rm} option set to TRUE.
}
\usage{
naMean(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An R object}
}
\value{
  The arithmetic mean of the values in \code{x}.
}
\examples{
  mean(c(1,2,3,NA),na.rm=TRUE)
  naMean(c(1,2,3,NA))
}
\keyword{ manip }
