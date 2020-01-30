\name{getUniprotBatch}
\alias{getUniprotBatch}
\title{
Function to query Uniprot database
}
\description{
Take a list of identifiers, assemble the Uniprot query, run and process the output
}
\usage{
getUniprotBatch(filters = "accession", values = IDList, attributes = c("id", "protein+names", "go-id"))
}
\arguments{
  \item{filters}{
 The Uniprot filter, typically the set of protein accessions
}
  \item{values}{
 The vector of protein ID's
}
  \item{attributes}{
 The fields to be extracted from Uniprot, typically protein id's, names and GO information
}
}
\details{
}
\value{
Data frame containing all extracted proteins, and their GO information.
}
\references{

}
\author{
D. Pascovici
}
\note{

}

\seealso{
genWegoFile
}
\examples{

values = c("Q9HWC9","Q9HWD0","Q9I4N8","Q9HW18","Q9HWC9","Q9HWD0")
getUniprotBatch(filters= "accession", values, attributes = c("id", "go-id"))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~utilities }
