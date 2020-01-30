\name{GOParent}
\alias{GOParent}
\title{
Function to find the first parent of a GO node
}
\description{
Find all parents and then select the first one.
}
\usage{
GOParent(node)
}
\arguments{
  \item{node}{
	GO node
}
}
%\details{
%}
\value{
Single GO node for the parent
}
%\references{
%}
\author{
D. Pascovici
}
%\note{
%}

\seealso{
\code{\link{Children}}
}
\examples{

GOParent("GO:0009636")
GOParent("bla")

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
