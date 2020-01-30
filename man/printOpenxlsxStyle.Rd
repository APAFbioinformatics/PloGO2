\name{printOpenxlsxStyle}
\alias{printOpenxlsxStyle}
\title{
Function for printing an Excel spreadsheet with highlighted cells
}
\description{
Function for printing an Excel spreadsheet with highlighted cells
for the p-values and ratios exceeding certain specified cut-offs.
}
\usage{
	printOpenxlsxStyle(dat, ratios, pvals, wb, tabName = "results", hiCutoff = 1.5, lowCutoff=0.67, pvalCutoff=0.05) 
}
\arguments{
\item{dat}{
  A data frame for print
}
\item{ratios}{
  Column indices identifying ratios
}
\item{pvals}{
  Column indices identifying pvalues
}
\item{wb}{
  openxlsx workbook object 
}
\item{tabName}{
   A name for the tab to be created
}
\item{hiCutoff}{
  Cutoff of upregulation FC
}
\item{lowCutoff}{
  Cutoff of downregulation FC
}
\item{pvalCutoff}{
  Cutoff of p value
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
\seealso{\code{\link{printSimpleOpenxlsxStyle}}
}
\references{
}
\examples{

wb <- openxlsx::createWorkbook("Example.xlsx")
data(iris)
printOpenxlsxStyle(iris[1:10,], ratios=c(1:2), pvals=4, hiCutoff=4, lowCutoff=3.1, pvalCutoff=0.3, tabName="Iris", wb = wb)
openxlsx::saveWorkbook(wb, file="Example.xlsx",  overwrite=TRUE)

}
\keyword{datasets}
