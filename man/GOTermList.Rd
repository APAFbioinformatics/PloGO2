\name{GOTermList}
\alias{GOTermList}
\title{
Function to extract a list of GO terms at level 2,3 or 4 of the GO hierarchy
}
\description{
At the moment a quick and dirty way to extract all GO nodes at levels 2,
3 or 4 of the GO hierarchy.  
}
\usage{
GOTermList(ontology = "BP", level = 2, node = NULL)
}
\arguments{
  \item{ontology}{
 Either "BP", "CC" or "MF"
}
  \item{level}{
 Either 2, 3 or 4
}
  \item{node}{
 Null or a GO node
}
}
\details{
This is a quick and dirty way to extract a list of GO terms of interest.
That can be one of the levels 2, 3 or 4 or all subnodes (children) of 
a particular node if a valid GO is provided for the node parameter.
Should be rewritten.
}
\value{
A vector of GO nodes
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

GOTermList("BP", 2)
GOTermList("CC", 2)


}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
