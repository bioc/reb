\name{fromRevIsh}
\alias{fromRevIsh}
\title{ Convert from revish strings to a matrix  }
\description{
  This function will convert two lists of revish style strings to a matrix format.
}
\usage{
fromRevIsh(enhList, dimList, chr, organism = "h")
}

\arguments{
  \item{enhList}{ list of enhanced bands on each individual sample }
  \item{dimList}{ list of diminished bands on each individual sample }
  \item{chr}{ chromosome to examine }
\item{organism}{character, "h" for human, "m" for mouse, and "r" for rat.}
}

\value{

A matrix is returned. The rownames of this matrix correspond to the major bands located on that chromosome, 
and the columns correspond to the sample names.

}
\references{ MCR eset data was obtained with permission. See PMID: 15377468 }
\author{ 
  Karl J. Dykema, \email{karl.dykema@vai.org} 
  Kyle A. Furge, \email{kyle.furge@vai.org}
}



\seealso{\code{\link{reb}},\code{\link{revish}} }
\examples{


mb.chr <- buildChromCytoband("h")

data(mcr.eset)
data(idiogramExample)
	## Create a vector with the index of normal samples
norms <- grep("MNC",colnames(mcr.eset@exprs))
	## Smooth the data using the default 'movbin' method, with the normal samples as reference and median centering
cset <- reb(mcr.eset,vai.chr,ref=norms,center=TRUE)
	## Mask the cset to remove noise
exprs <- cset@exprs[,-norms]
exprs[abs(exprs) < 1.96] <- NA
	## Extract the aberrations on the 5th chromosome
revish <- revish(exprs,vai.chr,"5")
	## Convert back to matrix
reconverted <- fromRevIsh(revish[[1]],revish[[2]],"5")


layout(cbind(1,2))
idiogram(cset@exprs[,-norms],vai.chr,"5",method="i",dlim=c(-2,2),col=.rwb,main="chr 5 reb results")
idiogram(reconverted,mb.chr,"5",method="i",dlim=c(-1,1),col=.rwb,main="chr 5 converted \n and re-converted")


	
}
\keyword{ manip}
