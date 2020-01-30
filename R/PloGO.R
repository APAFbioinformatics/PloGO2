
PloGO <- function(zipFile = "none", termFile="none", ontology = "BP", 
		ontologyLevel = 2, reference = "none",
		data.file.name = "none", datafile.ignore.cols = 1, 
		filesPath=".", node=NULL, aggregateFun="sum", logAb=FALSE, ...) 
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
	files <- unzip(zipFile)
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

if (DATA) {
	Group <- NULL
	gp <- names(read.csv(data.file.name))[-c(1:datafile.ignore.cols)]
	gp <- inferGroup(gp)
	if (min(summary(gp)) == max(summary(gp))) {
		Group <- gp;
		results$inferredGroup <- Group
	}
	abundance <- abundancePlot(res.list, Group=Group, log=logAb);
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


################
# Past examples
################

# PloGO()

# PloGO(termFile = "Z:/Projects/AnalysesCollection/Doc/sampleGOList/GOlistDrought.txt",
#	data.file.name = "Z:/Projects/AnalysesCollection/PloGO/inst/files/NSAFDesc.csv",
#	datafile.ignore.cols = 2  )


# PloGO(termFile = "Z:/Projects/AnalysesCollection/Doc/sampleGOList/GOlistDrought.txt",
#	data.file.name = "Z:/Projects/ExploreData/Steve-Grape/Grape8Conditions/Quality/NSAF.csv",
#	datafile.ignore.cols = 2, reference="Control"  )




# PloGO(data.file.name = "Z:/Projects/ExploreData/Prasanth/for merging dec11/Analysis IDS/Quality/NSAFDesc.csv",
#	datafile.ignore.cols = 2, reference="LungAll"  )





# PloGO(termFile = "Z:/Projects/AnalysesCollection/Doc/sampleGOList/GOlistDrought.txt",
#	data.file.name = "none", # "Z:/Projects/ExploreData/Iniga/Lorne poster/Quality/NSAF.csv",
#	datafile.ignore.cols = 1, reference="Control"  )





# PloGO(termFile = "Z:/Projects/AnalysesCollection/Doc/sampleGOList/GOlistDefense.txt",
#	data.file.name = "Z:/Projects/ExploreData/Steve-Grape/PloGO grant app/NSAFGrape.csv", 
#	datafile.ignore.cols = 2, reference="Ctr04"  )



# PloGO(termFile = "Z:/Projects/AnalysesCollection/Doc/sampleGOList/GOlistDefense.txt",
#	data.file.name = "Z:/Projects/AnalysesCollection/PloGO/inst/files/NSAFDesc.csv", 
#	datafile.ignore.cols = 2, reference="Control"  )


# PloGO(termFile = "Z:/Projects/AnalysesCollection/Doc/sampleGOList/GOlistDefense.txt",
#	data.file.name = "Z:/Projects/AnalysesCollection/PloGO/inst/files/NSAFDesc.csv", 
#	datafile.ignore.cols = 2, reference="Control", zipFile="fileZip.zip" )


# PloGO(termFile = "Z:/Projects/AnalysesCollection/Doc/sampleGOList/GOlistDefense.txt",
#	data.file.name = "Z:/Projects/AnalysesCollection/PloGO/inst/files/NSAFDesc.csv", 
#	datafile.ignore.cols = 2, reference="Control", zipFile="fileZip.zip" )


# res.list <-  PloGO(termFile = "Z:/Projects/AnalysesCollection/Doc/sampleGOList/GOGrapeExp.txt",
#	data.file.name = "Z:/Projects/ExploreData/Iniga/Lorne poster/Quality/NSAF.csv", 
#	datafile.ignore.cols = 2, reference="Control")


# PloGO(data.file.name = "Z:/Projects/ExploreData/Iniga/Lorne poster/Quality/NSAF.csv",
#	ontology="MF", ontologyLevel=3, 
#	datafile.ignore.cols = 2, reference="Control")



# PloGO( zipFile="fileZip.zip", reference="Control" )


