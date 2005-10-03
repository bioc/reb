\name{summarizeByRegion}
\alias{summarizeByRegion}
\alias{cgma}
\title{Compute Summary Statistics of Genome Regions}
\description{
Splits the data into subsets based on genome mapping information,
computes summary statistics for each region, and returns the results in
a convenient form. (cgma stands for Computer Generated Microarray Analysis)

This function supplies a t.test function at the emperically derived significance
threshold (p.value = 0.005)
}

\usage{
cgma(eset, genome, chrom="ALL",ref=NULL,center=TRUE,aggrfun=NULL, p.value=0.005, FUN=t.test, explode=FALSE ,...) 
}
\arguments{
  \item{eset}{an exprSet object}
  \item{genome}{an chromLocation object, such as on produced by buildChromLocation or buildChromMap}
  \item{chrom}{ a character vector specifying the chromosomes to analyze}
  \item{ref}{ a vector containing the index of reference samples from
    which to make comparisons. Defaults to NULL (internally referenced samples)} 
   \item{center}{ boolean - re-center gene expression matrix columns. Helpful if \code{ref} is used}
   \item{aggrfun}{ a function to summarizes/aggregates gene
      expression values that map to the same locations. If NULL, all values are
      included. Also see \code{absMax}}
  \item{p.value}{p.value cutoff or NA}
  \item{FUN}{ function by which to summarize the data}
  \item{explode}{boolean - explode summary matrix into a full expression set?}
  \item{\dots}{further arguments pass to or used by the function}
}
\details{
Gene expression values are separated into subsets that based on the
'chromLocation' object argument.  For example, buildChromMap can be used
to produce a 'chromLocation' object composed of the genes that populate
human chromosome 1p and chromosome 1q. The gene expression values from
each of these regions are extracted from the 'exprSet' and a summary
statistic is computed for each region.

\code{cgma} is most straightforwardly used to identify
regional gene expression biases when comparing a test sample to a
reference sample. For example, a number of simple tests can be used to
determine if a genomic region contains a disproportionate number of
positive or negative log transformed gene expression ratios. The
presence of such a regional expression bias can indicates an underlaying
genomic abnormality.

If multiple clones map to the same genomic locus the aggregate.by.loc
argument can be used to include a summary value for the overlapping
expression values rather then include all of the individual gene
expression values. For example, if 50 copies of the actin gene are on
a particular array and actin changes expression under a given condition, it may
appear as though a regional expression bias exists as 50 values in a
small region change expression.

\code{regmap} is usually the best way to plot results of this function. \code{ideogram} 
can also be used if you set the "explode" argument to TRUE.

\code{\link{buildChromLocation.2}} can be used to create a chromLocation object in 
which the genes can be divided a number of different ways. Separating the data by chromosome 
arm was the original intent. If you use \code{\link{buildChromLocation.2}} 
with the "arms" argument to build your chromLocation object, set the "chrom" argument 
to "arms" in this function. 



}
\value{
  \item{m}{A matrix of summary statistics}
}
\references{Crawley and Furge, Genome Biol. 2002;3(12):RESEARCH0075. Epub 2002 Nov 25. }
\author{Kyle A. Furge}
\note{}

\seealso{\code{\link{buildChromMap}},\code{\link{t.binom.test}},\code{\link{regmap}},\code{\link[ideogram]{buildChromLocation.2}}}
\examples{

## 
## NOTE: This requires an annotation package to work.
##       In this example packages "hu6800" and "golubEsets" are used.
##       They can be downloaded from http://www.bioconductor.org
##       "hu6800" is under MetaData, "golubEsets" is under Experimental Data.

if(require(hu6800) && require(golubEsets)) {
   cloc <- buildChromMap("hu6800",c("1p","1q","2p","2q","3p","3q"))

   ## For one-color expression data
   ## compare the ALL samples to the AML samples
   ## not particularly informative in this example

   aml.ix <- which(golubTrain$"ALL.AML" == "AML")
   bias <- cgma(eset=golubTrain,ref=aml.ix,genome=cloc)
   regmap(bias,col=.rwb) 
} else print("This example requires the hu6800 and golubEsets data
   packages.")

## A more interesting example

## The mcr.eset is a two-color gene expression exprSet
## where cytogenetically complex (MCR), 
## cytogenetically simple (CN) leukemia samples
## and normal control (MNC) samples were profiled against
## a pooled-cell line reference
## The MCR eset data was obtained with permission. See PMID: 15377468

## Notice the dimished expression on chromosome 5 in the MCR samples
## and the enhanced expression on chromosome 11
## This reflects chromosome gains and losses as validated by CGH

   data("mcr.eset")
   data(ideogramExample)
   norms <- grep("MNC",colnames(mcr.eset@exprs))
   bias <- cgma(mcr.eset,vai.chr,ref=norms)
   regmap(bias,col=topo.colors(50)) 
}
\keyword{manip}
