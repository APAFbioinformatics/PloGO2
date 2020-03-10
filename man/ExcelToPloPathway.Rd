\name{ExcelToPloPathway}
\alias{ExcelToPloPathway}
\title{
Function to add KEGG pathway annotation to an Excel spread sheet
}
\description{
Function to add KEGG pathway annotation to an Excel spread sheet
}
\usage{
ExcelToPloPathway(fname, colName = "Uniprot", compareWithReference = "none", DB.name = "pathwayDB.csv", data.file.name = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{
	Name of the Excel spread sheet to be annotated
}
  \item{colName}{
	The name of the column containing the protein identifiers
}
  \item{compareWithReference}{
	The name of the tab that serves as basis for enrichment comparison.
}
  \item{DB.name}{
	The DB file name for the pathway
}
  \item{data.file.name}{
	Abundance data for PloGO if any
}
}
\details{

}
\value{
A list object, with the following values:
\item{Counts}{
	The pathway counts matrix summarized for all files
}

\item{Percentages}{
	The pathway percentages matrix summarized for all files
}

\item{Abundance}{
	The pathway percentages matrix summarized for all files, ONLY generated if an
	abundance file was provided
}
\item{aggAbundance}{
	Aggregated abundance matrix by combination of tab and abundance file columns for all pathways
}
\item{FisherPval}{
	The Fisher p-values matrix summarized for all files, ONLY generated if
	a reference file was provided
}

\item{res.list}{
	The full list result of processAnnotation function
}
}
\references{

}
\author{

}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{

path <- system.file("files", package = "PloGO2")

res <- ExcelToPloPathway(file.path(path, "ResultsWGCNA_Input4PloGO2.xlsx"), 
	colName="Uniprot", compareWithReference="AllData", DB.name=file.path(path, "pathwayDB.csv"),
	data.file.name = file.path(path, "Abundance_data.csv") )

}

