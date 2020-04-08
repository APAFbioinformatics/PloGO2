
printSimpleOpenxlsxStyle  <- function (dat, tabName = "results" , wb ) 
{

	addWorksheet(wb, sheetName=tabName)
	altLines <- createStyle(fgFill = "lightblue")	
	writeData(wb, tabName, dat, keepNA=FALSE)
	
	for ( i in seq_len(ncol(dat))) addStyle(wb, tabName, style=altLines, rows = seq(1,nrow(dat)+1, 2), cols = i, gridExpand=TRUE);

	
   
}


