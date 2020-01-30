

PloPathway <- function(zipFile = "none",
		reference = "none",
		data.file.name = "none", 
		datafile.ignore.cols = 1, 
		filesPath=".", aggregateFun="sum", logAb=FALSE, ...) 
{

options(stringsAsFactors = FALSE)


results <- list()

if ( !(zipFile %in% c("", "none") ) ) {
	if ( !( "./PWFiles" %in% list.files()) ) dir.create("PWFiles")
	files <- unzip(zipFile, exdir="PWFiles")
	cat(files)
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


  # All unique pathway IDs
  
  allAnnotID = unlist(lapply(file.list, function(x) {v=read.table(x, sep="\t")[,2];
      sapply(v, function(y) strsplit(y, split="\\s")[[1]])}))
  
  AnnotIDlist <- names(sort(table(allAnnotID), decreasing = TRUE))


# get the protein IDs from each file
list.prot.ids = lapply(file.list, function(x) read.table(x, sep="\t")[,1])
names(list.prot.ids) = basename(file.list)
results$list.prot.ids = list.prot.ids


res.list <- processAnnotation(file.list, AnnotIDlist, data.file.name, 
	printFiles=FALSE, datafile.ignore.cols=datafile.ignore.cols, 
	aggregateFun=aggregateFun)

annot.res <- annotationPlot(res.list, plot=TRUE, trimzero=TRUE, type="pathway")
results$Counts <- annot.res$counts
results$Percentages <- annot.res$percentages

if(FALSE) {	
counts <- sapply(res.list, FUN = function(v) v$counts[, 1])
N <- sapply(res.list, FUN = function(v) v$N)
fnames <- sapply(res.list, FUN = function(v) v$fname)
colnames(counts) <- fnames
percentages <- 100 * sweep(counts, 2, N, FUN = "/")
pathTerms <- res.list[[1]]$counts[,2]
rownames(counts) <- pathTerms
rownames(percentages) <- pathTerms


results$Counts <- counts
results$Percentages <- percentages
}
Group <- NULL

if (DATA) {

	
	gp <- names(read.csv(data.file.name))[-c(1:datafile.ignore.cols)]
	

	#gp <- inferGroup(gp)
	Group <- gp
	
	if (min(summary(gp)) == max(summary(gp))) { # is this still needed to pathway?? 
		Group <- gp;
		results$inferredGroup <- Group
	}
	abundance <- abundancePlot(res.list, Group=Group, log=logAb, Plot=TRUE);
	results$Abundance <- abundance
}


if (useReference) {

comp <- compareAnnot(res.list, reference);
results$FisherPval <- comp

write.csv(comp, file="CompareFisher.csv")


}

# aggregated abundance

if (!is.null(Group)) {
  ag.list <- aggregateAbundance(res.list, Group=Group)
  fnames <- sapply(res.list, FUN = function(v) {
    v$fname
  })
  ag.mat <- ag.list[[1]]
  
  for (jj in 2:length(ag.list)) ag.mat <- rbind(ag.mat, ag.list[[jj]]);
  
  Files <- as.factor(rep(fnames, each=nlevels(as.factor(Group))))
  rownames(ag.mat) <- paste(rownames(ag.mat), Files)

  results$ag.list = ag.list
  results$aggAbundance = ag.mat
  results$Group = Group
  }
  results$res.list <- res.list



  results


}

