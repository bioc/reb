\name{smoothByRegion}
\alias{smoothByRegion}
\alias{reb}

\title{ reb }
\description{
  This function \dQuote{smooths} gene expression data to assist in the
  identification of regional expression biases. 
}
\usage{
reb(eset, genome, chrom = "ALL", ref = NULL, center = FALSE,
aggrfun=absMax, method = c("movbin", "supsmu", "lowess"), ...)
}
\arguments{
  \item{eset}{ the expression set to analyze }
  \item{genome}{ an associated chromLoc annotation object }
  \item{chrom}{ a character vector specifying the chromosomes to analyze}
  \item{ref}{ a vector containing the index of reference samples from
    which to make comparisons. Defaults to NULL (internally referenced samples}
  \item{center}{ boolean - re-center gene expression matrix columns. Helpful if
    \code{ref} is used}
    \item{aggrfun}{ a function to summarizes/aggregates gene
      expression values that map to the same locations. Defaults to the
      maximum absolute value \code{absMax}. If NULL, all values are
      included.}
  \item{method}{ smoothing function to use - either \code{"supmu"},
    \code{"lowess"} or \code{"movbin"}}
  \item{\dots}{ additional paramaters to pass along to the smoothing function }
}

\details{
  \code{reb} returns an eset that contains predictions of regional expression
  bias using data smoothing approachs. The exprSet is separated into
  subsets based on the \code{genome} chromLocation object and the gene
  expression data within the subsets is organized by genomic location
  and smoothed. In addition, the \code{approx} function is used to
  estimate data between any missing values. This was implimented so the
  function follows the \sQuote{principles of least astonishment}.    

  Smoothing approachs are most straightforwardly applied by comparing a
  set of test samples to a set of control samples.  For single color
  experiments, the control samples can be specified using the
  \code{ref} argument and the comparisons are generated internal to the
  \code{reb} function. This argument can also be used for two-color
  experiments provided both the test and control samples were run against
  a common reference.
  
  If multiple clones map to the same genomic locus the \code{aggrfun}
  argument can be used to summarize the overlapping expression
  values to a single summarized value. This is can be helpful in two
  situtations. First, the \code{supsum} and \code{lowess} smoothing
  functions do not allow for duplicate values. Currently, if duplicate
  values are found and these smoothing functions are used, the duplicate
  values are simply discard. Second, if 50 copies of
  the actin gene are present on a the array and actin changes expression
  under a given condition, it may appear as though a regional expression
  bias exists as 50 values within a region change expression.
  Summarizing the 50 expression values to a single value can partially
  correct for this effect. 
 
  The ideogram package can be used to plot the regional expression bias.
}

\value{
  An exprSet
}
\references{
Furge KA, Dykema KJ, Ho C, Chen X. Comparison of array-based comparative
genomic hybridization with gene expression-based regional expression
biases to identify genetic abnormalities in hepatocellular
carcinoma. BMC Genomics. 2005 May 9;6(1):67. PMID: 1588246

MCR eset data was obtained with permission. See PMID: 15377468
}
\author{ 
  Kyle A. Furge, \email{kyle.furge@vai.org}
  Karl J. Dykema, \email{karl.dykema@vai.org} 
}

\seealso{ \code{\link{movbin}},\code{ideogram}  }

\examples{

# The mcr.eset is a two-color gene expression exprSet
# with cytogenetically complex (MCR) and normal 
# control (MNC) samples which are a pooled-cell line reference.


data("mcr.eset")
data(ideogramExample)

## Create a vector with the index of normal samples
norms <- grep("MNC",colnames(mcr.eset@exprs))

## Smooth the data using the default 'movbin' method,
## with the normal samples as reference

cset <- reb(mcr.eset,vai.chr,ref=norms,center=TRUE)

## Display the results with mideogram
mideogram(cset@exprs[,-norms],vai.chr,method="i",dlim=c(-5,5),col=.rwb)
}
\keyword{manip}

