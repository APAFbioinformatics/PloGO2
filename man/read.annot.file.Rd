\name{read.annot.file}
\alias{read.annot.file}
\title{
Function to read an annotation file.
}
\description{
Accepts GO id's separated by space or semicolon
}
\usage{
read.annot.file(fname, format = c("compact","long"))
}
\arguments{
  \item{fname}{
	The file name containing GO annotation
}
  \item{format}{
	Either \code{"compact"} or \code{"long"}; see details
}
}
\details{
The format is \code{"compact"} by default, meaning a text file containing
ID's in the first column, and GO identifiers in the second, separated by
spaces or semicolons. This is the same as the "Wego native format".
A \code{"long"} format is also accepted, meaning a text file
with two or more columns separated by tabs, containing an identifier,
followed by a GO id, followed optionally by other columns which are ignored.
The GO id's will first be aggregated for each identifier.  
GO files in long format can be obtained using for instance biomart, 
or GoRetriever.

 
}
\value{
	A data frame with two columns, ID's and GO separated by spaces
}
%\references{
%}
\author{
T. Keighley, D.Pascovici
}
%\note{
%}


%\seealso{
%}
\examples{


# use one of the stored files
dir <- system.file("files", package="PloGO2")
fname <- paste(dir,"00100.txt", sep="/")

# Example with GoRetriever download
longFormat <- paste(dir,"goRetOutput.txt", sep="/")
read.annot.file(fname)
read.annot.file(longFormat, format="long")

# Example with biomart download
biomartDownload <- paste(dir,"mart_export.txt", sep="/")
read.annot.file(biomartDownload, format="long")



}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
