\name{extractPairs}
\alias{extractPairs}
\title{
Extract annotation and abundance side by side for two files 
}
\description{
Extract a data frame with the annotation and abundance for two chosen files.
}
\usage{
extractPairs(res.list, ag.list, pair = names(ag.list)[1:2], targetFile)
}
\arguments{
  \item{res.list}{
	Result of processAnnotation or processPlainAnnotation
}
  \item{ag.list}{
	Result of aggregateAbundance
}
  \item{pair}{
	Two of the file names used to generate the above
}
  \item{targetFile}{	
	One of the groups used when aggregating abundance
}
}
\details{
}
\value{
	A data frame with four columns containing annotation and abundance side by side
}
\references{
}
\author{
D.Pascovici
}
\note{
}

\seealso{
}
\examples{

# get list of ID's
termList <- c("response to stimulus", "transport", "protein folding","protein metabolic process", "carbohydrate metabolic process","oxidation reduction", "photosynthesis", "lipid metabolic process","cell redox homeostasis", "cellular amino acid and derivative metabolic process", "nucleobase, nucleoside and nucleotide metabolic process", "vitamin metabolic process", "generation of precursor metabolites and energy", "metabolic process", "signaling")
GOIDmap <- getGoID(termList)
GOIDlist <- names(GOIDmap)

# find existing files
dir <- system.file("files", package="PloGO2")
file.names <- paste(dir,c("00100.txt", "01111.txt", "10000.txt",
			"11111.txt","Control.txt"), sep="/")
datafile <- paste(dir, "/NSAF.csv", sep="")


# summarize annotation
res.list <- processAnnotation(file.names, GOIDlist, data.file.name=datafile)


# aggregate it by levels of the group
ag.list <- aggregateAbundance(res.list, Group=rep(c("Control", "Moderate", "Extreme", "3DaysRec", "6DaysRec"), each=3))


extractPairs(res.list, ag.list, c("00100.txt", "10000.txt"), "Extreme")



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
