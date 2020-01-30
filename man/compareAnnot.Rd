\name{compareAnnot}
\alias{compareAnnot}
\title{
Function to compare annotation percentages
}
\description{
Compare annotation percentages by means of Fisher's exact test.
A reference must be selected, and that name must be amongst the
annotation result list names.
}
\usage{
compareAnnot(res.list, referenceName, removeZeros = FALSE, correction = FALSE)
}
\arguments{
  \item{res.list}{
   List returned by the processAnnotation function
}
  \item{referenceName}{
   Name of the condition to compare with 
}
  \item{removeZeros}{
   Remove the categories with no annotation in them from the result
}
  \item{correction}{
   TRUE or FALSE: apply the BH fdr correction to the p-values in each column
}
}
%\details{
  
%}
\value{
 A matrix of p-values or NA.  A comparison is made for each sample other than the reference
 with the reference, and each GO annotation category.  The result is not recorded if the
 annotation numbers were small for that category (<5).  This is strictly not needed for
 Fisher's exact test (though was needed for the chi-square approximation used initially).

}
%\references{
%}
\author{
D.Pascovici
}
%\note{
%}


\seealso{
See Also \code{\link{processAnnotation}}
}
\examples{
\dontrun{
# get list of ID's
GOIDlist <- GOTermList("BP", 2)

# find existing files
dir <- system.file("files", package="PloGO2")
file.names <- paste(dir,c("00100.txt", "01111.txt", "10000.txt",
			"11111.txt","Control.txt"), sep="/")

# summarize annotation
res.list <- processAnnotation(file.names, GOIDlist)

# compare with Control
compareAnnot(res.list, "Control")

compareAnnot(res.list, "Control", correction=TRUE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
