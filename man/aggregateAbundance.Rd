\name{aggregateAbundance}
\alias{aggregateAbundance}
\title{
Function to aggregate the abundance by the group provided
}
\description{
Aggregates the rows of the abundance matrix by rows, giving one average value for each 
level of the grouping factor.  This can only be done if a data matrix was provided
when processing annotation. 
}
\usage{
aggregateAbundance(res.list, Group)
}
\arguments{
  \item{res.list}{
	The process annotation result 
}
  \item{Group}{
	A grouping factor of the same length as the number of data columns in the datafile provided to processAnnotation.
	
}
}
\details{

}
\value{
A list of the same length as the original res.list.  Each element is a matrix, containing averaged
data values for each level of the group.  Typically the group shows the technical replicates in an experiment,
so the data values will be averaged across the technical replicates. 
}
\references{
}
\author{
D. Pascovici
}
\note{
}


\seealso{
See Also \code{\link{processAnnotation}}
}
\examples{


# get list of ID's
GOIDlist <- GOTermList("MF", 2)

# find existing files
dir <- system.file("files", package="PloGO2")
file.names <- paste(dir,c("00100.txt", "01111.txt", "10000.txt",
			"11111.txt","Control.txt"), sep="/")
datafile <- paste(dir, "/NSAF.csv", sep="")


# summarize annotation
res.list <- processAnnotation(file.names, GOIDlist, data.file.name=datafile)

# aggregate it by levels of the group
ag.list <- aggregateAbundance(res.list, Group=rep(c("Control", "Moderate", "Extreme", "Recovery3", "Recovery6"), each=3))

}

