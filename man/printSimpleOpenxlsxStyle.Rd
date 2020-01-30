\name{printSimpleOpenxlsxStyle}
\alias{printSimpleOpenxlsxStyle}
\title{
Function for printing spreadsheet with alternating colours style
}
\description{
Utility function for printing simple alternating colour tables.
The spreadsheet must be created and closed before.
}
\usage{
printSimpleOpenxlsxStyle(dat, tabName = "results" , wb ) 
}
\arguments{
\item{dat}{
  A data frame for printing
}
\item{tabName}{
  A name for the tab to be created
}
\item{wb}{
  openxlsx workbook object 
}
}
\details{
This function is needed only for the printing of Excel summaries at the end.
Called from ExcelToPloGO.
}
\value{
 None returned
}
\source{
}
\note{
}
\seealso{\code{\link{printOpenxlsxStyle}}
}
\references{
}
\examples{


wb <- openxlsx::createWorkbook("Example.xlsx")
data(iris)
printSimpleOpenxlsxStyle(iris[1:10,], "Iris", wb = wb)
openxlsx::saveWorkbook(wb, file="Example.xlsx",  overwrite=TRUE)



}
\keyword{datasets}
