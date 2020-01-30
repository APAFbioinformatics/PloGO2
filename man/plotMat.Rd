\name{plotMat}
\alias{plotMat}
\title{
A levelplot of abundance by GO category
}
\description{
This plots a levelplot of the provided data matrix, sorted
by the number of proteins in each category 
}
\usage{
plotMat(z, number, cor = FALSE, log = FALSE, main = "Abundance Levelplot", ...)
}
\arguments{
  \item{z}{
	Data matrix
}
  \item{number}{
	A vector of numbers of length = nrow(z)
}
  \item{cor}{
	Boolean: plot correlation or data?
}
  \item{log}{
	Boolean: log before plotting?
}
  \item{main}{
	Plot title
}
\item{...}{
	Parameters to pass 
}
}
%\details{
%}
\value{
TRUE
}
%\references{
%}
\author{
D. Pascovici
}
%\note{
%}


%\seealso{
%}
\examples{
library(lattice)
z <- matrix(rnorm(100), nrow=20)

counts <- sample(1:100, 5)
plotMat(z, counts)

plotMat(z, counts, cor=TRUE)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
