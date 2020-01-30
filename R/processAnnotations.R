

processAnnotation <- function(file.list, AnnotIDlist, data.file.name = NULL, 
	printFiles = FALSE, format="compact", datafile.ignore.cols=1,
	aggregateFun="sum") 
  
  {
	annotation.list <- lapply(file.list, FUN = function(fname) {
	  
	  if(length(grep("^GO:", AnnotIDlist)) == length(AnnotIDlist)) {
	    res <- processGoFile(fname, AnnotIDlist, data.file.name, format=format, 
	                         datafile.ignore.cols=datafile.ignore.cols,
	                         aggregateFun=aggregateFun)
  	
	  } else {
	    res <- processPathFile(fname, AnnotIDlist, data.file.name, format=format, 
	                           datafile.ignore.cols=datafile.ignore.cols,
	                           aggregateFun=aggregateFun)
	  }

	  fout <- basename(fname)

		if (printFiles)
			writeGOannot(res, fname = paste("Annot", fout), datafile = data.file.name)

		res
	})

	names(annotation.list) <- basename(file.list)

	annotation.list
}
