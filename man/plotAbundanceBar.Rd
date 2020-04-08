\name{plotAbundanceBar}
\alias{plotAbundanceBar}
\title{
Function to plot abundance barplot
}
\description{
Plot the aggregated abundance barplot for all files/clusters

}
\usage{
plotAbundanceBar(mat.abundance, mat.counts, min.count=5)
}
\arguments{
\item{mat.abundance}{
  A matrix of abundance
}
\item{mat.counts}{
  A matrix of the counts
}
\item{min.count}{
  The cutoff for the minimum counts to be included 
}

}
%\details{
%}
\value{
 None returned, generate a barplot.
}
%\references{
%}
\author{
J.Wu
}
%\note{
%}


\examples{

path <- system.file("files", package = "PloGO2")
# run PloPathway with list of id's, data file and reference
res <- PloPathway( zipFile=paste(path, "PWFiles.zip", sep="/"), 
	reference="Alldata", 
	data.file.name = paste(path, "Abundance_data.csv", sep="/"),
	datafile.ignore.cols = 1)
	
png("AbundanceBarplot.png", 2500, 2000, res=300)
par(mar=c(4,10,4,14))

plot.res <- plotAbundanceBar(res$aggregatedAbundance, res$Counts)

dev.off()


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
