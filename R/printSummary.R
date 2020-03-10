
# New function - print summary and enrichment results

printSummary <- function(results, file="PloGO2Results.xlsx") {
  Counts = results$Counts
  colnames(Counts) = paste("Counts", colnames(Counts))
  Percentages = results$Percentages
  colnames(Percentages) = paste("Percentages", colnames(Percentages))
  
  res.summary = data.frame(Category=rownames(Counts), Counts, Percentages)
  if ( !is.null(results$FisherPval) ){
    FisherPval = results$FisherPval[rownames(results$FisherPval) %in% rownames(results$Counts),]
    colnames(FisherPval) = paste("FisherPval", colnames(FisherPval))
    res.summary = data.frame(res.summary, FisherPval)
  }
  
 
  
  annotated.list <- list()
  
  # make sure the folder order of files from PloGO is the
  # same as the sheet order from the spreadsheet!
  # sheet.ord <- match(paste(sheets, ".txt", sep=""), names(plogo.list[[1]]$res.list))
  
  # JW
  #sheet.ord <- match( names(results$res.list), paste(sheets, ".txt", sep=""))
  #data.list <- data.list[sheet.ord]
  
  data.list <- results$list.prot.ids
  
  
  for ( iii in 1:length(data.list) ) {
    
    prot.id <- data.list[[iii]]

    dat.res <- data.frame(Protein=prot.id, Order=1:length(prot.id))

    
    # JW
   adj.mat <- tabulateAnnot( results$res.list[[iii]] )
   
   dat.res <- merge( dat.res, adj.mat, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
   annotated.list[[iii]] <- dat.res[order(dat.res$Order),]

  }
  
  # Print results to tabbed sheet

  names(annotated.list) <- names(data.list)
  
  wb <- openxlsx::createWorkbook("Results.xlsx")
  
  
  for (jj in 1:length(annotated.list)) {
    if ("StouffersPval" %in% names(annotated.list[[jj]])) {
      pvals <- c(grep("Stouffer", names(annotated.list[[jj]])), 
                 grep("\\.pval\\.1..\\.1..$", tolower(colnames(annotated.list[[jj]]))))
      ratios <- c(grep("Geomean", names(annotated.list[[jj]])), 
                  grep("\\.X1..\\.1..$", colnames(annotated.list[[jj]])))
      ratios <- setdiff(ratios, grep("^Crit", colnames(annotated.list[[jj]])))
      tabName <- names(annotated.list)[jj]
    }
    else {
      ratios <- grep("^X1..\\.1..$", colnames(annotated.list[[jj]]))
      pvals <- grep("^pval\\.1..\\.1..$", tolower(colnames(annotated.list[[jj]])))
      tabName <- names(annotated.list)[jj]
    }
    printOpenxlsxStyle(annotated.list[[jj]], ratios, pvals, wb = wb, 
                       tabName = tabName)
  }
  
  
  tabName = "Summary"
  
  printSimpleOpenxlsxStyle(res.summary, "Summary", wb = wb)
  
  
  openxlsx::saveWorkbook(wb, file=file,  overwrite=TRUE)
  
}