\name{processPathFile}
\alias{processPathFile}
\title{
Function to process a single pathway file in wego-like native format
}
\description{
For each pathway extract all ID's from the file.
If abundance data is present extract and merge.
}
\usage{
processPathFile(fname, AnnotIDlist, datafile=NULL, datafile.ignore.cols=1, 
	format="compact", aggregateFun="sum")	
}
\arguments{
  \item{fname}{
	The pathway file name, in either Wego native format or long format
}
  \item{AnnotIDlist}{
	The list of pathway annotation ID
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
	The aggregation function for abundance data
}


}
\details{
The format is "compact" by default, meaning a text file containing
ID's in the first column, and pathway identifiers in the second, separated by
spaces or semicolons. This is the same as the "Wego native format".
A "long" format is also accepted, meaning a text file
with two or more columns separated by spaces, containing an identifier,
followed by a pathway id, followed optionally by other columns which are ignored.
The pathway id's will first be aggregated for each identifier.  
}
\value{
A list with the following components
\item{counts }{A vector of the same length as the list of pathway id's of interest giving the number of ID's in each category}
\item{ID.list}{The list of ID's for each pathway category }
\item{datafile}{The abundance datafile provided passed through}
\item{abundance}{A matrix with as many rows as the pathway list provided, and as many columns as the abundance data file }
\item{N}{The number of protein (gene etc) identifiers in each file}
\item{fname}{The filename without the file path}

}
%\references{
%}
\author{
J. Wu
}
%\note{
%}

%\seealso{
%}

\examples{

# use one of the stored files
dir <- system.file("files", package="PloGO2")
fname <- paste(dir,"PWFiles/Alldata.txt", sep="/")
datafile <- paste(dir, "Abundance_data.csv", sep="/")
AnnotIDlist <- unique(unlist(sapply(read.delim(fname, stringsAsFactors=FALSE)[,2], function(x) strsplit(x, split=" ")[[1]])))


# or if abundance in present aggregate that by category
processPathFile(fname, AnnotIDlist, datafile=datafile)


}


\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
