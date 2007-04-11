\name{movt}
\alias{movt}
\title{movt}
\description{
  This function analyzes ordered data series to identify regional
  biases using an moving (running) approximated t-test.
}
\usage{
movt(v,span=NULL,summarize=mean)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{data vector}
  \item{span}{numeric vector. Each element is used to define
    the number of points to include when the approximated binomial test
    is applied to \code{v}. While mixed for the defaults, the span can
    be specified as fraction of the observation or actual sizes, but
    \emph{not} a mixture - defaults to: seq(25,length(v)*.3,by=5)}
  \item{summarize}{function that is used to summarize the results from
    multiple spans. if NULL, a matrix with \code{length(span)} rows
    and \code{length(v)} columns is returned.}
}

\details{
  \code{movt} acts very similar to \code{movbin}
}

\value{
  Either a matrix or a vector containing the summarized z-scores from the applied t-test.
}

\author{ Kyle A. Furge, Ph.D., \email{kyle.furge@vai.org} and Karl J. Dykema, \email{karl.dykema@vai.org}}
\seealso{
\code{\link{movbin}}
}

\examples{
x <- c(rnorm(50,mean=1),rnorm(50,mean=-1),rnorm(100))
layout(1:2)
plot(x,type="h",ylim=c(-5,5))

## apply the approximated binomial with a single span
mb <- movbin(x,span=25,summarize=NULL)
lines(mb[1,])

## try a few different span ranges
mb <- movt(x,span=c(10,25,50),summarize=NULL)
lines(mb[1,])  ## span of 10
lines(mb[2,])  ## span of 25
lines(mb[3,])  ## span of 50

## average the results from the different spans
plot(x,type="h",ylim=c(-5,5))

mb <- movt(x,span=c(10,25,50),summarize=mean)
lines(mb,col="blue")

mb <- movt(x,span=c(10,25,50),summarize=median)
lines(mb,col="red")  

mb <- movt(x,span=c(10,25,50),summarize=max)
lines(mb,col="green")  

}
\keyword{manip}

