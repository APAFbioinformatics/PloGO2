\name{PloPathway}
\alias{PloPathway}
\title{
Function to do steps of pathway annotation and annotation/abundance plots 
}
\description{
Summarize pathway categories for all the files in the zip,
if provided merge data from data file,
generate annotation and abundance plots and comparison with reference.
}
\usage{
PloPathway(zipFile = "none", reference = "none", data.file.name = "none", datafile.ignore.cols = 1, filesPath = ".", aggregateFun="sum", logAb=FALSE,...)
}
\arguments{
  \item{zipFile}{
	Zip containing all the GO files
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
Process all the pathway files provided, in the directory or the zip,
and assign ID's to the respective categories.  
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
J. Wu
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{processAnnotation}}
}
\examples{
\dontrun{
# where sample files are stored
path <- system.file("files", package = "PloGO2")
# run PloPathway with list of id's, data file and reference
res <- PloPathway( zipFile=paste(path, "PWFiles.zip", sep="/"), 
	reference="Alldata", 
	data.file.name = paste(path, "Abundance_data.csv", sep="/"),
	datafile.ignore.cols = 1)
}
}
\keyword{ analysis }
