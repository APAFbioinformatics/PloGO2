\name{tabulateAnnot}
\alias{tabulateAnnot}
\title{
Function to create incidence matrix from list of GO categories
}
\description{
Utility to go from a list format of identifiers by ontology category to an incidence matrix of protein IDs as rows and
category as columns.
}
\usage{
tabulateAnnot(res.list)
}
\arguments{
  \item{res.list}{	
	Element of processAnnotation list
}
}
\details{
Utility called from ExcelToPloGO for summarising data prior to printing.
}
\value{
An incidence matrix with protein IDs as rows, and GO categories as columns.
}
\references{
}
\author{
D.Pascovici
}
\note{
}
\seealso{\code{\link{processGoFile}}
}
\examples{

termList <- c("response to stimulus", "transport")
GOIDmap <- getGoID(termList)
GOIDlist <- names(GOIDmap)

path <- system.file("files", package = "PloGO2")
file.names <- file.path(path, c("Control.txt", "00100.txt"))

# summarize annotation
res.list <- processAnnotation(file.names, GOIDlist)

tabulateAnnot(res.list[[1]])

}
\keyword{ manip }
