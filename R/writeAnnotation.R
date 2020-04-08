writeAnnotation <- function(res.list, datafile=NULL, datafile.ignore.cols=1, format=c("list","matrix"), outFolder=tempdir()) {

format <- match.arg(format)
res <- lapply(res.list, FUN=function(res) {
	fout <- basename(res$fname)
	if (format == "matrix") fout <- gsub("txt", "csv", fout);
	writeGOannot(res, fname = (paste("Annot", fout)), datafile = datafile, format=format, outFolder=outFolder)

})

outFolder

}


