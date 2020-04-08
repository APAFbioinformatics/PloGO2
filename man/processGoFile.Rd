\name{processGoFile}
\alias{processGoFile}
\title{
Function to process a single file in wego-like native format
}
\description{
For each GO category of interest extract all ID's from the file
that are annotated at either the category or its GO children.
If abundance data is present extract and merge.
}
\usage{
processGoFile(fname, GOIDlist, datafile = NULL, datafile.ignore.cols = 1, format = c("compact","long"),aggregateFun="sum")
}
\arguments{
  \item{fname}{
	The GO file name, in either Wego native format or long format
}
  \item{GOIDlist}{
	The list of GO id's of interest
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
 \item{aggregateFun}{
	Either "sum" or "product"; the aggregation operation for abundance data
}
}
\details{
The format is "compact" by default, meaning a text file containing
ID's in the first column, and GO identifiers in the second, separated by
spaces or semicolons. This is the same as the "Wego native format".
A "long" format is also accepted, meaning a text file
with two or more columns separated by spaces, containing an identifier,
followed by a GO id, followed optionally by other columns which are ignored.
The GO id's will first be aggregated for each identifier.  The output of
GOretriever can be used as "long" format.
}
\value{
A list with the following components
\item{counts }{A vector of the same length as the list of GO id's of interest giving the number of ID's in each category}
\item{ID.list}{The list of ID's for each GO category }
\item{datafile}{The abundance datafile provided passed through}
\item{abundance}{A matrix with as many rows as the GO list provided, and as many columns as the abundance data file }
\item{N}{The number of protein (gene etc) identifiers in each file}
\item{fname}{The filename without the file path}

}
%\references{
%}
\author{
D. Pascovici
}
%\note{
%}

%\seealso{
%}

\examples{

termList <- c("response to stimulus", "transport")
GOIDmap <- getGoID(termList)
GOIDlist <- names(GOIDmap)

# use one of the stored files
dir <- system.file("files", package="PloGO2")
fname <- paste(dir,"00100.txt", sep="/")
datafile <- paste(dir, "NSAF.csv", sep="/")


# or if abundance in present aggregate that by category
processGoFile(fname, GOIDlist, datafile=datafile)



}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
