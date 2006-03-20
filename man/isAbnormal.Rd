\name{isAbnormal}
\alias{isAbnormal}
\title{ Is a band 'abnormal'?}
\description{
	Returns 1 or -1 indicating a chromosomal change based upon an input percentage.
}
\usage{
isAbnormal(x, percent = 0.5)
}
\arguments{
  \item{x}{ genomic data, can contain NA's }
  \item{percent}{ numeric argument - a fraction or percentage}
}
\details{
  This simple function is used by \code{cset2band}. 
}
\value{
  
}

\author{ Karl Dykema  }

\seealso{ \code{\link{cset2band}} }
\examples{
		#Not abnormal
	isAbnormal(c(1,NA))
		#Abnormal; +
	isAbnormal(c(1,NA,1))
		#Abnormal; -
	isAbnormal(c(1,NA,-1,-1,-1))
	

}
\keyword{ arith }

