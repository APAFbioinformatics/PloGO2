\name{processPlainAnnotation}
\alias{processPlainAnnotation}
\title{
Function to process a set of annotation files given a list of identifiers of interest
}
\description{
For each file in the list extract all identifiers that belong to each
category in the list of identifiers.  
}
\usage{
processPlainAnnotation(files, IDlist, datafile = NULL, printFiles = FALSE, format = "compact", datafile.ignore.cols = 1, term.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
	A list of files
}
  \item{IDlist}{
	A list of identifiers
}
  \item{datafile}{
	A list containing additional data such as abundance information
}
  \item{printFiles}{
	TRUE/FALSE If true an annotation summary file is generated for each input file
}
  \item{format}{
	One of "compact" by default or or "long", the format of the annotation files.  See details.
}
  \item{datafile.ignore.cols}{
	How many columns in the abundance file to ignore.  By default assume the first only, containing identifiers.
}
  \item{term.names}{
	A data frame with two columns, identifiers and their names
}
}
\details{
The format is for the annotation files is "compact" by default, meaning a text file containing
ID's in the first column, and identifiers in the second, separated by
spaces or semicolons. This is the same as the "Wego native format".
A "long" format is also accepted, meaning a text file
with two or more columns separated by spaces, containing an identifier,
followed by an id, followed optionally by other columns which are ignored.
The id's will first be aggregated for each identifier.  
}
\value{
A list of the same length as the list of files submitted.  Each list element is
	the result of processAnnotFile.
}
\references{
}
\author{
D. Pascovici
}
\note{
}

\seealso{
\code{\link{processAnnotFile}} \code{\link{processGoFile}} 
}
\examples{

# choose two simple GO categories
IDlist <- c("J", "D")
term.names <- data.frame(IDlist, 
	c("Translation, ribosomal structure and biogenesis", 
		"Cell cycle control, cell division, chromosome partitioning"))


dir <- system.file("files", package="PloGO2")
file.names <- paste(dir,c("test1.txt", "test2.txt"), sep="/")
datafile <- paste(dir,"NSAFtest.csv", sep="/")


# summarize annotation
res.list <- processPlainAnnotation(file.names, IDlist)

res.list <- processPlainAnnotation(file.names, IDlist, datafile=datafile, term.names=term.names)

# or infer the list 
res.list <- processPlainAnnotation(file.names, datafile=datafile)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
