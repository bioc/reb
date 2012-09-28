\name{rmAmbigMappings}
\alias{rmAmbigMappings}

\title{Remove genes that map to multiple chromosomes from a chromLocation object}
\description{
  Due to the automated probe annotation, a subset of probes can be
  \dQuote{confidently} mapped to multiple chromosomes on the genome.

  This can cause some confusion if you are trying to perform certain
  types of data analysis.

  This function examines a \code{chromLocation} object and removes
  probes that map to multiple chromosomes.
}
\usage{
rmAmbigMappings(cL)
}
\arguments{
  \item{cL}{an existing \code{chromLocation} object
}
}
\value{
  A \code{chromLocation} object
}

\author{Kyle A. Furge}

\seealso{\code{\link[annotate]{buildChromLocation}}}

\examples{

## 
## NOTE: This requires an annotation package to work, it is provided for info only.
## 

#if (require(hu6800)) {

#  library(Biobase)
#  library(annotate)

  ## Build a specific chrom arm

#  cL <- buildChromLocation("hu6800")
#  cleanCL <- rmAmbigMappings(cL)
#}
}

\keyword{ manip }
