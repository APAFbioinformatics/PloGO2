\name{abundancePlot}
\alias{abundancePlot}
\title{
Function to summarize and plot abundance information from an annotations results list
}
\description{
Generates one GO/pathway abundance plot for each file, and (provided
the number of GO/pathway cateGO/pathwayries is reasonably small) also an abundance
plot for each GO/pathway category across all files provided.
The abundance values are *added* for all values in one category.
Hence if the initial values represented percentages such as NSAF,
the final values represent percentages of the respective category.

}
\usage{
abundancePlot(res.list, log = FALSE, printLimit = 16, Group=NULL, Plot=FALSE, 
	CountCutOff=3, ...)
}
\arguments{
\item{res.list}{
  The result returned by processAnnotation
}
\item{log}{
  TRUE/FALSE: use raw or log abundance data
}
\item{printLimit}{
  The most number of category to be plotted
}
\item{Group}{
  The groups
}
\item{Plot}{
  To plot or not
}
\item{CountCutOff}{
  The minimum number of proteins in a category for it to be plotted
}
\item{...}{
	Parameters to pass
}
}
%\details{
%}
\value{
A list object, with the following values:
\item{abundance}{
	The abundance matrix
}
\item{ag.mat}{
	The aggregated abundance matrix
}
\item{list.levelplots}{
	A list of abundance levelplots trellis object
}

\item{list.barplots}{
	A list of abundance barcharts trellis object
}
 

}
%\references{
%}
\author{
D.Pascovici and J.Wu
}
%\note{
%}


\seealso{
See Also as \code{\link{processAnnotation}}
}
\examples{

# get list of ID's
GOIDlist <- GOTermList("MF", 2)

# find existing files
path <- system.file("files", package="PloGO2")
file.names <- paste(path, c("00100.txt", "01111.txt", "10000.txt",
			"11111.txt","Control.txt"), sep="/")
datafile <- file.path(path, "NSAF.csv")

# summarize annotation
res.list <- processAnnotation(file.names, GOIDlist, data.file.name = datafile)

abundance.res <- abundancePlot(res.list)

# Plot levelplots
list.levelplots <- abundance.res$list.levelplots

for(i in seq_along(list.levelplots)) {
	png(paste(names(list.levelplots)[i], ".png"), 2000, 4000, res=200)
	print(list.levelplots[[i]])
	dev.off()
}	


# KEGG pathway
path <- system.file("files", package="PloGO2")
file.names <- file.path(path,"PWFiles", c("AllData.txt","black.txt","blue.txt","brown.txt","green.txt",
	"red.txt","turquoise.txt") )
datafile <- file.path(path,"Abundance_data.csv")
Group <- names(read.csv(datafile))[-1]

AnnotIDlist <- c("osa01100","osa01110","osa01230","osa00300","osa00860")
		
res.list <- processAnnotation(file.names, AnnotIDlist,  data.file.name = datafile)

abundance.res <- abundancePlot(res.list, Group=Group, Plot=TRUE)

# Plot levelplots
list.levelplots <- abundance.res$list.levelplots

for(i in seq_along(list.levelplots)) {
	png(paste(names(list.levelplots)[i], ".png"), 2000, 4000, res=200)
	print(list.levelplots[[i]])
	dev.off()
}	
# Plot barchats by categories
list.barplots <- abundance.res$list.barplots

for(i in seq_along(list.barplots)) {
	png(paste(names(list.barplots)[i], ".png"), 2000, 2000, res=200)
	print(list.barplots[[i]])
	dev.off()
	
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
