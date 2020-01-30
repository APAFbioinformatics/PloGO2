\name{inferGroup}
\alias{inferGroup}
\title{
Utility to infer grouping by cleaning up vector of names, typically deriving from file names
}
\description{
Remove dots and ending csv filenames as well as replicate number
}
\usage{
inferGroup(v)
}
\arguments{
  \item{v}{
	Vector of file names
}
}
\details{
Utility only
}
\value{
Cleaned up vector of names
}
\references{
}
\author{
}
\note{
}
\seealso{
}
\examples{
v = c("Control.1.csv", "Control.2.csv", "Control.3.csv", "Stress.1.csv", "Stress.2.csv", "Stress.3.csv")
inferGroup(v)
}
\keyword{misc}
