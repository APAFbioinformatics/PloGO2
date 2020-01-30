\name{countAndAbundance}
\alias{countAndAbundance}
\title{
Function to do a side by side annotation and abundance plot
}
\description{
Plots annotation and abundance values side by side.
}
\usage{
countAndAbundance(x, main = "Annotation and abundance plot", legend.text = NULL, args.legend = list(x = "topright"))
}
\arguments{
  \item{x}{
	A matrix with annotation numbers in the first half, and abundance in the second half
}
  \item{main}{
	The image title
}
  \item{legend.text}{
	The legend text
}
  \item{args.legend}{
	Arguments for the legend
}
}
\details{
}
\value{
The plot generated
}
\references{
}
\author{
T.Keighley
}
\note{
}

\seealso{
\code{\link{extractPairs}}
}
\examples{

# get list of ID's
termList <- c("response to stimulus", "transport", "protein folding","protein metabolic process", "carbohydrate metabolic process","oxidation reduction", "photosynthesis", "lipid metabolic process","cell redox homeostasis", "cellular amino acid and derivative metabolic process", "nucleobase, nucleoside and nucleotide metabolic process", "vitamin metabolic process", "generation of precursor metabolites and energy", "metabolic process", "signaling")
GOIDmap <- getGoID(termList)
GOIDlist <- names(GOIDmap)

# find existing files
dir <- system.file("files", package="PloGO2")
file.names <- paste(dir, c("00100.txt", "01111.txt", "10000.txt",
			"11111.txt","Control.txt"), sep="/")
datafile <- paste(dir,"/NSAF.csv", sep="")


# summarize annotation
res.list <- processAnnotation(file.names, GOIDlist, data.file.name=datafile)


# aggregate it by levels of the group
ag.list <- aggregateAbundance(res.list, Group=rep(c("Control", "Moderate", "Extreme", "3DaysRec", "6DaysRec"), each=3))


x <- extractPairs(res.list, ag.list, c("00100.txt", "10000.txt"), "Extreme")

countAndAbundance(x, "Two random files", legend.tex=c("Extreme Only", "Control Only"))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
