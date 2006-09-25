\name{cset2band}
\alias{cset2band}

\title{ cset2band}
\description{
  This function will summarize gene expression data by cytogenetic band
}
\usage{
cset2band(exprs, genome, chr = "ALL", organism = NULL, FUN = isAbnormal, ...)
}

\arguments{
  \item{exprs}{ matrix of gene expression data or similar. The rownames
  must contain the gene identifiers }
  \item{genome}{ an associated chromLoc annotation object }
  \item{chr}{a character vector specifying the chromosomes to analyze}
  \item{organism}{ character, "h" for human, "m" for mouse, and "r" for rat.; defaults to NULL - loads from chromLocation object }
  \item{FUN}{ function by which to aggregate/summarize each cytogenetic band }
  \item{\dots}{ extra arguments passed on to the aggregate/summary function }
}
\details{
  This function loops through each band for a given organism and
  summarizes the data for genes that lie  within each cytogenetic band
  based upon the input function.
  For example, a matrix of gene expression values could be used
  and the mean expression of each band be determined by passing the
  \code{mean} function.
  Alternative, DNA copy number gains or losses could be predicted using
  the \code{reb} function and regions of likely gain or losses be
  summarized by cytogenetic band using the \code{isAbnormal} function.
}
\value{
  a matrix with rows representing cytogenetic bands, and columns representing individual samples.
}

\author{ Karl Dykema}


\examples{

   data(mcr.eset)
   data(idiogramExample)

## Create a vector with the index of normal samples
   norms <- grep("MNC",colnames(mcr.eset@exprs))

## Smooth the data using the default 'movbin' method,
## with the normal samples as reference and median centering
   cset <- reb(mcr.eset,vai.chr,ref=norms,center=TRUE)

## Mask the result to remove noise
   exprs <- cset[,-norms]
   exprs[abs(exprs) < 1.96] <- NA

## Starting data
   midiogram(exprs,vai.chr,method="i",col=.rwb,dlim=c(-4,4))

## Summarize each cytogenetic band
   banded <- cset2band(exprs,vai.chr,FUN=mean,na.rm=TRUE)

## Create chromLocation object based on human cytobands
   h.cyto <- buildChromCytoband(organism = "h")

## Plot all data using mideogram
   midiogram(banded,h.cyto,method="i",col=.rwb,dlim=c(-4,4))

}
\keyword{ manip }

