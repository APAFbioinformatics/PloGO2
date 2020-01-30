\name{Children}
\alias{Children}
\title{
Function to extract the children 
}
\description{
Extracts all children of a GO node irrespective of type of relationship "is_a" etc. Utility for GOTermList.
}
\usage{
Children(node)
}
\arguments{
  \item{node}{
	A string representing a valid GO node
}
}
%\details{
%}
\value{
An array of GO nodes, or NULL if the node is not a proper GO node
}
%\references{
%}
\author{
D. Pascovici
}
%\note{
%}
\seealso{
See Also \code{\link{GOTermList}}
}
\examples{

Children("GO:0009636")
Children("bla") 

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
