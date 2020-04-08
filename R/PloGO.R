
PloGO <- function(zipFile = "none", termFile="none", ontology = "BP", 
		ontologyLevel = 2, reference = "none",
		data.file.name = "none", datafile.ignore.cols = 1, 
		filesPath=".", node=NULL, aggregateFun="sum", logAb=FALSE,...) 
{


results <- list()

if (!(termFile %in% c("", "none"))) {
	termList <- as.vector(read.delim(termFile)[,1])
	tpList <- getGoID(termList)
	GOIDlist <- names(tpList)
	results$termList <- termList
} else {
	GOIDlist <- GOTermList(ontology, as.numeric(ontologyLevel), node)
	results$GOLevel <- ontologyLevel
	results$GOOntology <- ontology
}


if ( !(zipFile %in% c("", "none") ) ) {
  # unzip to tempdir
	files <- unzip(zipFile, exdir=tempdir())
	#cat(files)
} else {
	files <- list.files(filesPath)
 	files <- paste(filesPath, files, sep="/")
}



file.list <- files[grep("\\.txt$", files)]
# don't try to analyze Annot files just in case from previous run
if (length(-grep("^Annot", file.list)) > 0) file.list <- file.list[-grep("^Annot", file.list)];
file.names <- gsub("\\.txt", "", unlist(lapply(file.list, FUN=basename)))


DATA <- TRUE
if (data.file.name %in% c("", "none") ) {
	DATA <- FALSE
	data.file.name <- NULL
} 

results$datafileName <- data.file.name
results$datafile.ignore.cols <- datafile.ignore.cols

useReference <- (!(reference == "none")) 


# JW, get the protein IDs from each file
list.prot.ids = lapply(file.list, function(x) read.annot.file(x)[,1])
names(list.prot.ids) = basename(file.list)
results$list.prot.ids = list.prot.ids


res.list <- processAnnotation(file.list, GOIDlist, data.file.name, 
	printFiles=FALSE, datafile.ignore.cols=datafile.ignore.cols, 
	aggregateFun=aggregateFun)

# should review annotationPlot - DP Nov 2019
annot.res <- annotationPlot(res.list, plot=TRUE, trimzero=TRUE)
results$Counts <- annot.res$counts
results$Percentages <- annot.res$percentages
Group <- NULL

if (DATA) {
	gp <- names(read.csv(data.file.name))[-seq_len(datafile.ignore.cols)]
	gp <- inferGroup(gp)
	if (min(summary(gp)) == max(summary(gp))) {
		Group <- gp;
		results$inferredGroup <- Group
	}
	abundance <- abundancePlot(res.list, Group=Group, log=logAb);
	results$Abundance <- abundance$abundance
	results$aggregatedAbundance <- abundance$ag.mat
	
	# The levelplots and barcharts
	results$list.levelplots <- abundance$list.levelplots
	results$list.barplots <- abundance$list.barplots
}


if (useReference) {

comp <- compareAnnot(res.list, reference);
results$FisherPval <- comp

}

results$Group = Group  
results$res.list <- res.list

results


}

