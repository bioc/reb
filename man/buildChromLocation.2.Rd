\name{buildChromLocation.2}
\alias{buildChromLocation.2}
\title{A function to generate an instantiation of a chromLocation class }
\description{
  This function will take the name of a data package and build a
  chromLocation object representing that data set. It has also been modified to allow further breakup of the chromLocs.
}
\usage{
  buildChromLocation.2(dataPkg,major=NULL)
}
\arguments{
  \item{dataPkg}{The name of the data package to be used}
  \item{major}{name of major breakpoint by which to divide chromosomes, "arms", "bands", and "mb" currently work.}
}
\details{
  The requested data set must be available in the user's
  \code{.libPaths()}, and the function will throw an error if this is
  not the case.

  If the data package is present, the necessary information will be
  extracted from the data package and a \code{chromLocation} object will
  be created.
  
  If "major" is set to "arms", the the chromLocs object is populated with data from the chromosome arms; "1p", "1q", "2p", etc... Rat and Human chromsomes follow this pattern , so data packages from both species should work with this function.

  If "major" is set to "bands", the chromosomes are divided up based upon which band they fall into. 

  If "major" is set to "mb", chromosomes are split into 3000+ megabase segments. Note, this creates a large chromLocation object.
  
  Note, "major' can contain multiple breakpoint names, eg. major=c("arms","bands")
  
}
\value{
  A \code{chromLocation} object representing the specified data set.
}

\author{ Jeff Gentry and Karl Dykema}

\examples{
   ## A bit of a hack to not have a package dependency on hgu95av2
   ## but need to fiddle w/ the warn level to not fail the example anyways.
   curWarn <- getOption("warn")
   options(warn=0)
   on.exit(options(warn=curWarn), add=TRUE)
   if (require(hgu95av2) & require(ideogram)) {
     z <- buildChromLocation.2("hgu95av2",major="arms")
   } else print("This example requires the hgu95av2 data package")
}
\keyword{ utilities }
