\name{plotAbundance}
\alias{plotAbundance}
\title{
Plot abundance data by category and file/cluster
}
\description{
Plot abundance data by category and file/cluster
}
\usage{
	plotAbundance(res.list, log=FALSE, eps = 0.001)
}
\arguments{
\item{res.list}{
  The result returned by processAnnotation
}
\item{log}{
  TRUE/FALSE: use raw or log abundance data
}
\item{eps}{
  A small value added to all if taking logs
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
J.Wu
}
%\note{
%}


\examples{
\dontrun{
file.names <- file.path(path,"PWFiles", c("black.txt","blue.txt","brown.txt","green.txt",
	"red.txt","turquoise.txt") )
datafile <- file.path(path,"Abundance_data.csv")
Group <- names(read.csv(datafile))[-1]

AnnotIDlist <- c("osa01100","osa01110","osa01230","osa00300","osa00860")
		
res.list <- processAnnotation(file.names, AnnotIDlist,  data.file.name = datafile)

plotAbundance(res.list)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
