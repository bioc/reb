\name{movbin}
\alias{movbin}
\title{movbin}
\description{
  This function analyzes ordered data series to identify regional
  biases using an moving (running) approximated binomial test.
}
\usage{
movbin(v,span=seq(25,length(v)*.3,by=5),summarize=mean)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{data vector}
  \item{span}{numeric vector. Each element is used to define
    the number of points to include when the approximated binomial test
    is applied to \code{v}. While mixed for the defaults, the span can
    be specified as fraction of the observation or actual sizes, but
    \emph{not} a mixture}
  \item{summarize}{function that is used to summarize the results from
    multiple spans. if NULL, a matrix with \code{length(span)} rows
    and \code{length(v)} columns is returned.}
}

\details{
  \code{movbin} applies a moving binomial test to sequential windows of
  elements of \code{v}. Within each span a z-score from an approximated
  binomial is computed such that \code{z=(2*r - n)/sqrt(n)} where
  \code{r} is the number of positive relative gene expression values and
  \code{n} is the number of non-zero values within each window.

  For convenience, this function allows for the specification of multiple
  window sizes using the \code{span} argument. The result of a
  \code{movbin} call will generate a matrix with \code{length(span)} rows
  and \code{length(v)} columns. Each row of the matrix represents the
  data generated from each span. This matrix can be returned or the matrix
  from can be condensed to a single vector of length v by
  applying a summary function \code{summarize} to the matrix columns.
}

\value{
  Either a matrix or a vector containing the summarized z-scores from the applied binomial test.
}

\author{ Kyle A. Furge, Ph.D., \email{kyle.furge@vai.org} and Karl J. Dykema, \email{karl.dykema@vai.org}}
\note{}

\examples{
x <- c(rnorm(50,mean=1),rnorm(50,mean=-1),rnorm(100))
layout(1:2)
plot(x,type="h",ylim=c(-5,5))

## apply the approximated binomial with a single span
mb <- movbin(x,span=25,summarize=NULL)
lines(mb[1,])

## try a few different span ranges
mb <- movbin(x,span=c(10,25,50),summarize=NULL)
lines(mb[1,])  ## span of 10
lines(mb[2,])  ## span of 25
lines(mb[3,])  ## span of 50

## average the results from the different spans
plot(x,type="h",ylim=c(-5,5))

mb <- movbin(x,span=c(10,25,50),summarize=mean)
lines(mb,col="blue")

mb <- movbin(x,span=c(10,25,50),summarize=median)
lines(mb,col="red")  
}
\keyword{manip}

