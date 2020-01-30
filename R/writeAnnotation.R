writeAnnotation <- function(res.list, datafile=NULL, datafile.ignore.cols=1, format="list") {


res <- lapply(res.list, FUN=function(res) {
	fout <- basename(res$fname)
	if (format == "matrix") fout <- gsub("txt", "csv", fout);
	writeGOannot(res, fname = paste("Annot", fout), datafile = datafile, format=format)

})


}


# writeAnnotation(res.list, format="matrix")
