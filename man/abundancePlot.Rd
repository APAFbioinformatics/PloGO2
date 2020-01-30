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
 None returned, generate several images in the local folder.
}
%\references{
%}
\author{
D.Pascovici
}
%\note{
%}


\seealso{
See Also as \code{\link{processAnnotation}}
}
\examples{
\dontrun{

# get list of ID's
GOIDlist <- GOTermList("MF", 2)

# find existing files
path <- system.file("files", package="PloGO2")
file.names <- paste(path, c("00100.txt", "01111.txt", "10000.txt",
			"11111.txt","Control.txt"), sep="/")
datafile <- file.path(path, "NSAF.csv")


# summarize annotation
res.list <- processAnnotation(file.names, GOIDlist, data.file.name = datafile)

ab.mat <- abundancePlot(res.list)

# abundance plot
image(ab.mat)
}


# KEGG pathway
path <- system.file("files", package="PloGO2")
file.names <- file.path(path,"PWFiles", c("Alldata.txt","black.txt","blue.txt","brown.txt","green.txt",
	"red.txt","turquoise.txt") )
datafile <- file.path(path,"Abundance_data.csv")
Group <- names(read.csv(datafile))[-1]

AnnotIDlist <- c("osa01100","osa01110","osa01230","osa00300","osa00860")
		
res.list <- processAnnotation(file.names, AnnotIDlist,  data.file.name = datafile)

ab.mat <- abundancePlot(res.list, Group=Group)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
