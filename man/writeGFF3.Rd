\name{writeGFF3}
\alias{writeGFF3}
\title{ Output of a GFF compliant table describing the enhanced and diminished chromosomal bands.}
\description{
  This function writes out a GFF compliant tab delimited file for intergration with genome browsers.
}
\usage{
writeGFF3(cset, genome, chr, file.prefix = "temp.gff", organism = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cset}{ expression set containing cytogenetic predictions, see \code{\link{reb}}  }
  \item{genome}{ chromLocation object containing annotation information }
  \item{chr}{ chromosome to examine }
  \item{file.prefix}{ character string - name of the output file, defaults to "temp.gff" }
\item{organism}{ if NULL, determination of the host organism will be
  retrieved from the \code{organism} slot of the chromLocation
  object. Otherwise \code{"h"}, \code{"r"}, or \code{"m"} can be used to
  specify \bold{h}uman, \bold{r}at, or \bold{m}ouse chromosome
  information}
}

\value{
\code{writeGFF3} returns an invisible list of character vectors.
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
ix <- abs(colo.eset@exprs) > .225
colo.eset@exprs[ix] <- NA
idiogram(colo.eset@exprs,ucsf.chr,"14",method="i",dlim=c(-1,1),col=.rwb)
gffmat <- writeGFF3(colo.eset@exprs,ucsf.chr,"14",NULL)
gffmat[1:4,]

}
\keyword{ manip}
