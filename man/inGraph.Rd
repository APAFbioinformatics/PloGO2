\name{inGraph}
\alias{inGraph}
\title{
Function to decide if a particular term is in a graph
}
\description{
Utility
}
\usage{
inGraph(terms, graph)
}
\arguments{
  \item{terms}{
	A GO term
}
  \item{graph}{
	A GO graph as returned by GOGraphWrapper
}
}
%\details{
%}
\value{
TRUE or FALSE
}
%\references{
%}
\author{
D. Pascovici
}
%\note{
%}

\seealso{
\code{\link{GOGraph}}
}
\examples{

g <- GOGraphWrapper("GO:0005576")
inGraph("GO:0005576", g)
inGraph("zz", g)

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
