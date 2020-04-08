\name{writeAnnotation}
\alias{writeAnnotation}
\title{
Function to print GO/pathway annotation to files
}
\description{
Prints available GO or pathway annotation and abundance (if existing) in
a long format or an adjacency matrix type format.
}
\usage{
writeAnnotation(res.list, datafile = NULL, datafile.ignore.cols = 1, format = c("list","matrix"),outFolder=tempdir())
}
\arguments{
  \item{res.list}{
	The result of processAnnotation
}
  \item{datafile}{
	A CSV file with additional experimental information, if any
}
  \item{datafile.ignore.cols}{
	The number of columns to ignore in the data file
}
  \item{format}{
  	Either \code{"matrix"} or \code{"list"}.
}
 \item{outFolder}{
	The output files folder
}
}
\details{
The GO/pathway information and abundance will be printed to files.  If the format
is \code{"list"}, then the files will be text files, and each category will
be printed in turn, with all the identifiers and data underneath.
If the format is \code{"matrix"}, then the data will be printed in matrix format,
identifiers (rows) by GO categories (columns), with the abundance data appended.
}
\value{
 The path of the annotation folder.
}
%\references{
%}
\author{
D. Pascovici
}
%\note{
%}
\seealso{
\code{\link{processAnnotation}}
}
\examples{

# choose two simple GO categories
termList <- c("response to stimulus", "transport", "signaling")
GOIDmap <- getGoID(termList)
GOIDlist <- names(GOIDmap)


dir <- system.file("files", package="PloGO2")
file.names <- paste(dir,c("00100.txt", "01111.txt", "10000.txt",
			"11111.txt","Control.txt"), sep="/")

# summarize annotation
res.list <- processAnnotation(file.names, GOIDlist)

# write to "matrix" or alternatively "list" format
writeAnnotation(res.list, format="matrix")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
