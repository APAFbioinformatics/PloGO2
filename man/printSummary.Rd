\name{printSummary}
\alias{printSummary}
\title{
Function to print the summary file of PloGO2 results.
}
\description{
Each tab in the summary spreadsheet of the  
}
\usage{
	printSummary(results, file="PloGO2Results.xlsx") 
}
\arguments{
  \item{results}{
	A list of results from PloGO2 analysis 
} 
\item{file}{
	The output file name
} 

}

%\details{
%}
\value{
TRUE if at least some annotations were found, FALSE otherwise.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
J. Wu
}
%\note{
%}

\examples{
\dontrun{
path <- system.file("files", package = "PloGO2")
# run PloPathway with list of id's, data file and reference
res <- PloPathway( zipFile=paste(path, "PWFiles.zip", sep="/"), 
	reference="Alldata", 
	data.file.name = paste(path, "Abundance_data.csv", sep="/"),
	datafile.ignore.cols = 1)

printSummary(res)
}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

