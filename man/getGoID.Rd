\name{getGoID}
\alias{getGoID}
\title{
Function to map a vector of GO terms to the corresponding GO ID's
}
\description{
This function is rather slow and inefficient, as it
first extracts all GO terms then matches the current one.
However, it in only intended for a one off matching
of terms of interest to the respective nodes.

}
\usage{
getGoID(v)
}
\arguments{
  \item{v}{ the precise spelling of the go term}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The list of matched terms, with names the respective GO ID's/
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

getGoID(c("biological_process", "transport"))

# however the next one is not found as the proper term name has an 
# underscore "cellular_component"
getGoID("cellular component")

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
