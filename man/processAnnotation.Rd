\name{processAnnotation}
\alias{processAnnotation}
\title{
Function to process a set of annotation files given a list of GO/pathway identifiers of interest 
}
\description{
For each file in the list extract all identifiers that belong to each
GO or pathway category in the list of identifiers. An identifier \dQuote{belongs} to
a GO or pathway category if it is annotated at the category itself or any
of its children in the GO graph.

}
\usage{
processAnnotation(file.list, AnnotIDlist, data.file.name = NULL, printFiles = FALSE,
	format = "compact", datafile.ignore.cols = 1, aggregateFun="sum")
}	
\arguments{
  \item{file.list}{
	A list of files
}
  \item{AnnotIDlist}{
 	A list of GO/pathway identifiers
}
  \item{data.file.name}{
	A list containing additional data such as abundance information
}
  \item{printFiles}{
	TRUE/FALSE If true an annotation summary file is generated for each input file
}
  \item{format}{
	One of \dQuote{compact} by default or or \dQuote{long}, the format of the
	annotation files. See details.
}
  \item{datafile.ignore.cols}{
	How many columns in the abundance file to ignore.  By default assume the first only, containing identifiers.
}
 \item{aggregateFun}{
	Either "sum" or "product"; the aggregation operation for abundance data
}
}
\details{
The format is for the GO/pathway files is \dQuote{compact} by default, meaning a text file containing
ID's in the first column, and GO/pathway identifiers in the second, separated by
spaces or semicolons. This is the same as the \dQuote{Wego native format}.
A \dQuote{long} format is also accepted, meaning a text file
with two or more columns separated by spaces, containing an identifier,
followed by a GO/pathway id, followed optionally by other columns which are ignored.
The GO/pathway id's will first be aggregated for each identifier. 
}
\value{
	A list of the same length as the list of files submitted.  Each list element is
	the result of \code{\link{processGoFile}} or \code{\link{processPathFile}}.
}
%\references{
%}
\author{
D. Pascovici, J.Wu
}
%\note{
%}

\seealso{
\code{\link{processGoFile}}, \code{\link{processPathFile}}
}
\examples{

# For GO analysis
# choose two simple GO categories
termList <- c("response to stimulus", "transport")
GOIDmap <- getGoID(termList)
GOIDlist <- names(GOIDmap)


path <- system.file("files", package = "PloGO2")
file.names <- file.path(path, c("00100.txt", "01111.txt", "10000.txt",
	"11111.txt", "Control.txt"))

# summarize annotation
res.list <- processAnnotation(file.names, GOIDlist)

# files in long format
longFormat <- file.path(path,  "goRetOutput.txt")
processAnnotation(longFormat, GOIDlist, format = "long")

# For KEGG pathway analysis
fname <- file.path(path,"PWFiles", "red.txt")
datafile <- file.path(path,"Abundance_data.csv")

AnnotIDlist <- c("osa01100","osa01110","osa01230","osa00300","osa00860")
		
res.list <- processAnnotation(fname, AnnotIDlist,  data.file.name = datafile)


}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
