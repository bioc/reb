\name{revish}
\alias{revish}
\title{ Creation of CGH (reverse in situ hybridization) style character strings  }
\description{
	This function returns a two lists of character strings. These two lists correspond to 
the enhanced and diminished chromosomal bands.  
}
\usage{
revish(cset, genome, chr, organism = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cset}{ expression set containing cytogenetic predictions, see \code{\link{reb}}  }
  \item{genome}{ chromLocation object containing annotation information }
  \item{chr}{ chromosome to examine }
\item{organism}{ if NULL, determination of the host organism will be
  retrieved from the \code{organism} slot of the chromLocation
  object. Otherwise \code{"h"}, \code{"r"}, or \code{"m"} can be used to
  specify \bold{h}uman, \bold{r}at, or \bold{m}ouse chromosome
  information}
}


\value{
   
  \item{enh }{list of enhanced bands on each individual sample}
  \item{dim }{list of diminished bands on each individual sample}
  
}
\references{
Furge KA, Dykema KJ, Ho C, Chen X. Comparison of array-based comparative
genomic hybridization with gene expression-based regional expression
biases to identify genetic abnormalities in hepatocellular
carcinoma. BMC Genomics. 2005 May 9;6(1):67. PMID: 1588246

MCR eset data was obtained with permission. See PMID: 15377468
}
\author{ 
  Karl J. Dykema, \email{karl.dykema@vai.org} 
  Kyle A. Furge, \email{kyle.furge@vai.org}
}

\seealso{ \code{\link{reb}} }
\examples{

data(idiogramExample)
ix <- abs(colo.eset) > .225
colo.eset[ix] <- NA
idiogram(colo.eset,ucsf.chr,"14",method="i",dlim=c(-1,1),col=.rwb)
revlist<- revish(colo.eset,ucsf.chr,"14")
str(revlist)

}
\keyword{ manip}
