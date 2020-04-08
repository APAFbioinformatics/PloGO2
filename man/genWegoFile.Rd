\name{genWegoFile}
\alias{genWegoFile}
\title{
Function to generate a format in Wego native style from a list of Uniprot identifiers
}
\description{
Connects to the Uniprot Biomart using functionality from the biomaRt package,
downloads GO information and organizes it as needed.
}
\usage{
genWegoFile(IDList, fname = "Wego.txt", database = "uniprot", outFolder=tempdir())
}
\arguments{
  \item{IDList}{
	A list of ID's, either Ensembl or Uniprot.
}
  \item{fname}{
	The name of the text file to be outputed.
}
\item{database}{
	\code{"ensembl"} or \code{"uniprot"}.
}
 \item{outFolder}{
	The output files folder
}
}

%\details{
%}
\value{
File path of generated Wego files
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
D. Pascovici
}
%\note{
%}

\examples{
v <- c("Q9HWC9","Q9HWD0","Q9I4N8","Q9HW18","Q9HWC9","Q9HWD0")

\dontrun{genWegoFile(v, fname = "F1.txt")}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
