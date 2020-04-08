\name{ExcelToPloGO}
\alias{ExcelToPloGO}
\title{
Function to add GO annotation to an Excel spread sheet
}
\description{
Function to add GO annotation to an Excel spread sheet
}
\usage{
ExcelToPloGO(fname, colName = "Uniprot", termFile= NA, compareWithReference="none", data.file.name = "none", outFolder = "PloGO2Output")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{
	Name of the Excel spread sheet to be annotated
}
  \item{colName}{
	The name of the column containing the protein identifiers
}
  \item{termFile}{
	The name of the file containing the GO categories
}
  \item{compareWithReference}{
	The name of the tab that serves as basis for enrichement comparison.
}
  \item{data.file.name}{
	Abundance data for PloGO if any
}
 \item{outFolder}{
	The output files folder
}
}

\value{
A list object, with the following values:
\item{Counts}{
	The GO counts matrix summarized for all files
}

\item{Percentages}{
	The GO percentages matrix summarized for all files
}


\item{Abundance}{
	The GO percentages matrix summarized for all files, ONLY generated if an
	abundance file was provided
}

\item{FisherPval}{
	The Fisher p-values matrix summarized for all files, ONLY generated if
	a reference file was provided
}

\item{res.list}{
	The full list result of processAnnotation function
}



\item{list.levelplots}{
	The abundance levelplots if the data.file.name is not "none"
}

\item{list.barplots}{
	The abundance barcharts if the data.file.name is not "none"
}


}

\author{
D.Pascovici
}

\examples{
\dontrun{
# where sample files are stored
path <- system.file("files", package = "PloGO2")
termFile = paste(path, "GODefault.txt", sep="/")
xlfile <- paste(path, "ResultsWGCNA_Input4PloGO2.xlsx", sep="/")

res <- ExcelToPloGO(xlfile, termFile=termFile, compareWithReference="AllData")


}
}