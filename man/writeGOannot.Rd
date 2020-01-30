\name{writeGOannot}
\alias{writeGOannot}
\title{
Function to write annotation out to a file
}
\description{
Write one file summarizing the identifiers present
in each GO category for each sample (GO file).

}
\usage{
writeGOannot(res, fname = "AnnotOut.txt", datafile = NULL, format = "list")
}
\arguments{
  \item{res}{
   Annotation result list, produced by processWegoFile
}
  \item{fname}{
   Output file name
}
  \item{datafile}{
   Data file containing abundance results or NULL if none
}
  \item{format}{
  	Either \code{"matrix"} or \code{"list"}.
}
}
%\details{

%}
\value{
None returned, output file generated
}
%\references{
%}
\author{
D.Pascovici
}

%\note{
%}


\seealso{
See also \code{\link{processGoFile}}
}
\examples{

dir <- system.file("files", package="PloGO2")
fname <- paste(dir,"/00100.txt", sep="")
datafile <- paste(dir,"/NSAF.csv", sep="")

writeGOannot(processGoFile(fname, GOTermList("BP", 2)),
		datafile=datafile)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
