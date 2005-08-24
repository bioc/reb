\name{regmap}
\alias{regmap}
\title{image function wrapper}
\description{
A simple wrapper around the image function
}
\usage{
regmap(m,scale=c(-6,6),na.color=par("bg"),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{a matrix}
  \item{scale}{Include a graph scale showing this range of values \sQuote{image} function}
  \item{na.color}{the color to draw over NA values}
  \item{\dots}{additional paramters to \sQuote(image}}
}
\details{
A small wrapper around the \sQuote{image} function to display genome region summary statistics. Additional parameters will be passed along to the image function.  

The scale argument is a two-element vector that provides a floor and
ceiling for the matrix and allows a crude scale bar to be included on the lower
potion of the graph. 

For other colors consider using the geneplotter (dChip.colors) or marrayPlots (maPalette) library functions (i.e. regmap(m,col=dChipColors(50)))
}
\value{
}
\references{}
\author{Kyle A. Furge}
\note{}

\seealso{\code{\link{image}},\code{\link{summarizeByRegion}}}
\examples{
   m <- matrix(rnorm(6*4),ncol=6)
   colnames(m) <- c(1:6)
   rownames(m) <- c("1p","1q","2p","2q")
   regmap(m,scale=c(-1,1))
}
\keyword{manip}

