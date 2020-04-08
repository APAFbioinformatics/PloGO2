\name{PloGO}
\alias{PloGO}
\title{
Function to do steps of GO annotation and annotation/abundance plots 
}
\description{
Summarize GO categories for all the files in the zip,
if provided merge data from data file,
generate annotation and abundance plots and comparison with reference.
}
\usage{
PloGO(zipFile = "none", termFile="none", ontology = "BP", ontologyLevel = 2, reference = "none", data.file.name = "none", datafile.ignore.cols = 1, 
		filesPath=".", node=NULL, aggregateFun="sum", logAb=FALSE, ...) 
}
\arguments{
  \item{zipFile}{
	Zip containing all the GO files
}
  \item{termFile}{
	A file with the GO terms of interest, if a limited set provided
}
  \item{ontology}{
	Wither BP, MF or CC
}
  \item{ontologyLevel}{
	A small level such as 2 or 3
}
  \item{reference}{
	The file name of the reference file, if any, for instance "Control" for Control.txt
}
  \item{data.file.name}{
	The file containing all the experimental data, for instance "NASF.csv"
}
  \item{datafile.ignore.cols}{
	The number of columns in the experimental file given in data.file.name that should not be used as numerical,
		for instance 2 if there is an ID field and a Description field. By default 1.
}
  \item{filesPath}{
	If the zip file is not provided, the path to the GO files
}
  \item{node}{
	NULL, ignored at the moment
}
 \item{aggregateFun}{
	Either "sum" or "product"; the aggregation operation for abundance data
}
 \item{logAb}{
	TRUE or FALSE; the abundance data to be logged or not
}
 \item{...}{
	Parameters to pass
}
}
\details{
Process all the GO files provided, in the directory or the zip,
and assign ID's to the respective categories.  The categories can
come from a target list, or from a choice of level and ontology.
If a reference is provided, then the numbers of identifiers in
each category is compared to the reference by means of Fisher's
exact test.  If a data file is provided, then the quantitative
values are aggregated into the sets.

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

\item{aggregatedAbundance}{
	The aggregated abundance matrix if the data file is provided
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
D. Pascovici
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{processAnnotation}}
}
\examples{
# where sample files are stored
path <- system.file("files", package = "PloGO2")
# run PloGO with list of id's, data file and reference
res <- PloGO( zipFile=paste(path, "GOfiles.zip", sep="/"), 
	reference="Control", termFile = paste(path, "GOListDrought.txt", sep="/"),
	data.file.name = paste(path, "NSAFDesc.csv", sep="/"),
	datafile.ignore.cols = 2)

}
\keyword{ analysis }
