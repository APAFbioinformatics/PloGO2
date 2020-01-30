printOpenxlsxStyle  <- function (dat, ratios, pvals, wb, tabName = "results", hiCutoff = 1.5, lowCutoff=0.67, pvalCutoff=0.05) 
{

	addWorksheet(wb, sheet=tabName)
	
	upReg <- createStyle(fgFill = "violet")
	downReg <- createStyle(fgFill = "forestgreen")
	sigStyle <- createStyle(fgFill = "khaki1")
	
	writeData(wb, tabName, dat, keepNA=FALSE)
	
    for (rat in ratios) {
        up.idx <- which(!is.na(dat[, rat]) & (dat[, rat] > hiCutoff))
        if (length(up.idx) > 1) 
			addStyle(wb, tabName, style=upReg, rows = 1 + up.idx, cols = rat)

        down.idx <- which(!is.na(dat[, rat]) & (dat[, rat] < 
            lowCutoff))
        if (length(down.idx) > 1) 
			addStyle(wb, tabName, style=downReg, rows = 1 + down.idx, cols = rat)
    }
	
    for (pval in pvals) {
        sig.idx <- which(!is.na(dat[, pval]) & (dat[, pval] < 
            pvalCutoff))
        if (length(sig.idx) > 1) 
			addStyle(wb, tabName, style=sigStyle, rows = 1 + sig.idx, cols = pval)
    }
}
