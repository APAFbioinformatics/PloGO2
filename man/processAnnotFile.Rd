\name{processAnnotFile}
\alias{processAnnotFile}
\title{
Function to process a single file in wego-like native format
}
\description{
For each GO category of interest extract all ID's from the file
that are annotated at either the category or its GO children.
If abundance data is present extract and merge.
}
\usage{
processAnnotFile(fname, GOIDlist, datafile = NULL, datafile.ignore.cols = 1, format = "compact", term.names = NULL)
}
\arguments{
  \item{fname}{
	The GO file name, in either Wego native format or long format
}
  \item{GOIDlist}{
	The list of id's of interest
}
  \item{datafile}{
	The file containing abundance or NULL if none. 
}
  \item{datafile.ignore.cols}{
	How many columns in the abundance file to ignore.  By default assume the first only, containing identifiers.
}
  \item{format}{
	Either "compact" or "long"; see details
}
  \item{term.names}{
	A data frame with two columns, identifiers and their names
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The format is for the annotation files is "compact" by default, meaning a text file containing
ID's in the first column, and identifiers in the second, separated by
spaces or semicolons. This is the same as the "Wego native format".
A "long" format is also accepted, meaning a text file
with two or more columns separated by spaces, containing an identifier,
followed by an id, followed optionally by other columns which are ignored.
The id's will first be aggregated for each identifier.  
}
\references{
}
\author{
D. Pascovici
}
\note{
}

\seealso{
 \code{\link{processPlainAnnotation}},  \code{\link{processGoFile}}
}
\examples{

IDlist <- c("J", "D")
term.names <- data.frame(IDlist, 
	c("Translation, ribosomal structure and biogenesis", 
		"Cell cycle control, cell division, chromosome partitioning"))


dir <- system.file("files", package="PloGO2")
fname <- paste(dir,"test1.txt", sep="/")
datafile <- paste(dir,"NSAFtest.csv", sep="/")


# summarize annotation

res <- processAnnotFile(fname, IDlist, datafile=datafile, term.names=term.names)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
