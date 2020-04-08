# Author J. Wu, jemma.wu@mq.edu.au
# Date 28 Jan 2020
# Description:
# This script generates Figure 4 and Support information file 3 in 
# paper "Quickly characterise complex, multi-condition proteomics experiments with WGCNA and PloGO2"

	


library(PloGO2)

########################################
# ExcelToPloPathway for the main dataset
 ########################################
 
path <- system.file("files", package = "PloGO2")

res <- ExcelToPloPathway(file.path(path,"ResultsWGCNA_Input4PloGO2.xlsx"), 
	colName="Uniprot", compareWithReference="AllData", DB.name=file.path(path,"pathwayDB.csv"),
	data.file.name = file.path(path,"Abundance_data.csv") )



#############
# Paper plots	
#############

# Run the WGCNA_proteomics.R for Figure 2, WGCNA R package has to be installed first

# source(file.path(system.file("script", 	package="PloGO2"), "WGCNA_proteomics.R"))


# Level plots in Figure 4a and 4b

path <- system.file("files", package="PloGO2")
file.names <- file.path(path,"PWFiles", list.files(file.path(path,"PWFiles"), pattern=".txt") )
datafile <- file.path(path,"Abundance_data.csv")
Group <- names(read.csv(datafile))[-1]
options(stringsAsFactors=FALSE)
allAnnotID = unlist(lapply(file.names, function(x) {v=read.table(x, sep="\t")[,2];
  sapply(v, function(y) strsplit(y, split="\\s")[[1]])}))

AnnotIDlist <- names(sort(table(allAnnotID), decreasing = TRUE))
	
res.list <- processAnnotation(file.names, AnnotIDlist,  data.file.name = datafile)

list.levelplots <- 
	abundancePlot(res.list, Group=Group, CountCutOff=2, cex.main=2, cex.lab=1)$list.levelplots

for(i in seq_along(list.levelplots)) {
	png(paste("File",i,".png"), 2000, 4000, res=200)
	print(list.levelplots[[i]])
	dev.off()
	}
	
	

# Abundance barplot in Figure 4c

path <- system.file("files", package = "PloGO2")

datafile <- file.path(path,"Abundance_data.csv")


file.names <- file.path(path,"PWFiles", list.files(file.path(path,"PWFiles"), pattern=".txt") )
options(stringsAsFactors=FALSE)
allAnnotID = unlist(lapply(file.names, function(x) {v=read.table(x, sep="\t")[,2];
  sapply(v, function(y) strsplit(y, split="\\s")[[1]])}))

AnnotIDlist <- names(sort(table(allAnnotID), decreasing = TRUE))	
	
res.list <- processAnnotation(file.names, AnnotIDlist,  data.file.name = datafile)

res.annot = annotationPlot(res.list, plot=TRUE, type="pathway")

GOTerm <- rep(rownames(res.annot$percentages), ncol(res.annot$percentages))

Files <- rep(basename(file.names), each = nrow(res.annot$percentages))
idx=5
Percentage = res.annot$percentages[order(res.annot$percentages[,6], decreasing=TRUE),][1:idx,]

pdf("Annotation1Perc.pdf")
print(barchart(log(Percentage + 1) ~ rep(rownames(Percentage),ncol(Percentage)), 
	groups = rep(basename(file.names), idx),
			scales = list(x = list(rot = 60)), auto.key = list(points = FALSE,
			rectangles = TRUE, space = "top"), cex.axis=2))
dev.off()

# PloPathway for the main dataset
res <- PloPathway( zipFile=paste(path, "PWFiles.zip", sep="/"), 
reference="Alldata", 
data.file.name = paste(path, "Abundance_data.csv", sep="/"),
datafile.ignore.cols = 1)

idx.dataset = grep("blue|green|red|yellow|turquoise", rownames(res$aggAbundance))
idx.countcol = grep("blue|green|red|yellow|turquoise", colnames(res$Counts))
idx.pathway = grep(paste(rownames(Percentage), collapse="|"), colnames(res$aggAbundance) )

# Abundance barplot	
png("abundanceBarplot.png", 2500, 2000, res=300)
par(mar=c(4,10,4,14))
plot.res <- plotAbundanceBar(res$aggAbundance[idx.dataset, idx.pathway], 
	res$Counts[idx.pathway, idx.countcol], min.count=0)
dev.off()







