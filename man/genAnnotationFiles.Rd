\name{genAnnotationFiles}
\alias{genAnnotationFiles}
\title{
Function to generate a format in Wego native style from a list of Uniprot identifiers
}
\description{
Given a Excel spreadsheet with multiple tabs, generate a Wego file for each tab using a pre-download 
DB file.  
}
\usage{
genAnnotationFiles(fExcelName,  colName="Uniprot", 
				DB.name = "pathwayDB.csv", folder="PWFiles")
}
\arguments{
  \item{fExcelName}{
	An excel file containing one or multiple tabs of protein IDs. For example, the proteins in each tab
	come from the same cluster.
}
  \item{colName}{
	The column name of the protein ID.
}
\item{DB.name}{
	The database file name, in .csv format.
}
\item{folder}{
	The folder name for saving the generated files.
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

path <- system.file("files", package = "PloGO2")

genAnnotationFiles(fExcelName = file.path(path, "ResultsWGCNA_Input4PloGO2.xlsx"),  
				colName="Uniprot", 
				DB.name = file.path(path, "pathwayDB.csv"))

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
