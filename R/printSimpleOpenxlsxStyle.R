
printSimpleOpenxlsxStyle  <- function (dat, tabName = "results" , wb ) 
{

	addWorksheet(wb, sheet=tabName)
	altLines <- createStyle(fgFill = "lightblue")	
	writeData(wb, tabName, dat, keepNA=FALSE)
	
	for ( i in 1:ncol(dat)) addStyle(wb, tabName, style=altLines, rows = seq(1,nrow(dat)+1, 2), cols = i, gridExpand=TRUE);

	
   
}


