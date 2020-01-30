# Author J. Wu, jemma.wu@mq.edu.au
# Date 28 Jan 2020
# Description:
# This script generates Figure 2 and Support information file 1 in 
# paper "Quickly characterise complex, multi-condition proteomics experiments with WGCNA and PloGO2"

	
rm(list=ls())

library(openxlsx)
library(heatmap3)
library(WGCNA)

options(stringsAsFactors = FALSE)

# set working directory
# setwd("C:\\tmp")




############
# Functions
############

	

plotErrorBarsLines <- function (v, barSizes, lines, labels = NULL, col = "blue", 
	ylim = c(min(lines), max(lines)), ...) 
{
  barSizes[is.na(barSizes)] <- 0
  topBars <- v + 0.5 * barSizes
  bottomBars <- v - 0.5 * barSizes
  N <- length(v)
  if (is.null(labels)) 
    labels <- 1:N
  ylims <- c(min(bottomBars, ylim[1], min(lines)), max(topBars, 
                                                       ylim[2], max(lines)))
  par(pch = 19, xaxt = "n")
  plot(as.numeric(labels), v, ylim = ylims, col = col, type = "b", 
       lwd = 3, ...)
  par(xaxt = "s")
  
  for (i in 1:N) {
    lines(c(i, i), c(topBars[i], bottomBars[i]))
  }
  for (i in 1:ncol(lines)) {
    lines(as.numeric(labels), lines[, i], lwd = 0.5, lty = "dotted", 
          col = "gray")
  }
}

plotClusterProfileWGCNA <- function(cluster.data, moduleColors, group, MEs=NULL, 
		ylab="Abundance", 
		file="ClusterPatterns.png", ...) {
	
	gp = group
	noClusters <- nlevels(as.factor(moduleColors))
	
	r.temp <- aggregate(t(cluster.data), by=list(gp=gp), FUN=mean)
	ag.sample <- r.temp[,-1]
	rownames(ag.sample) <- r.temp[,1]
	ag.genes <- aggregate(t(ag.sample), by=list(Cluster=moduleColors), FUN=mean)
	ag.sd <- aggregate(t(ag.sample), by=list(Cluster=moduleColors), FUN=sd)
	ag.matrix <- as.matrix(ag.genes[,-1])

	if(!is.null(MEs) ) {		
		r.temp <- aggregate(MEs, by=list(gp=gp), FUN=mean)
		ag.matrix <- t(r.temp[,-1])
		colnames(ag.matrix) <- r.temp[,1]
	}
	ag.counts <- summary(as.factor(moduleColors))
	ag.bars <- as.matrix(ag.sd[,-1])
	
	fScale = max(8,noClusters)/8
	
	
	png(file, 2000, 3000*fScale, res=300)
	par(bg=gray(.95), fg=gray(0.3), mar= c(8, 6, 2, 1) + 0.1, col.main="black", col.sub="black", col.lab="black", col.axis="black")
	layout(matrix(1:(ceiling(noClusters/2)*2), ncol=2, byrow=TRUE))
	NSig <- noClusters
	cols = levels(as.factor(moduleColors) )
	for(i in 1:NSig) {
		gname <-  paste(levels(as.factor(moduleColors))[i], "(", ag.counts[i], "proteins )")
		lines <- ag.sample[, moduleColors==levels(as.factor(moduleColors))[i], drop=FALSE]
		plotErrorBarsLines(ag.matrix[i,], 2*ag.bars[i,], lines, 
			labels=1:ncol(ag.matrix), 
			col=cols[i],  main=gname, # bgcol="gray", split=split,
			ylab=ylab, xlab="",
			ylim=c(min(ag.matrix), max(ag.matrix)), ...)
		axis(1,at=1:ncol(ag.matrix), las=2, labels=colnames(ag.matrix), col="black", ...)
		abline(h=0, lty="dotted")
	}
	
	dev.off()

}


###############
# WGCNA workflow
###############

# Read the proteomics abundance data - it's assumed preprocessed and normalisation done

path <- system.file("files", package="PloGO2")
allDat = read.csv(file.path(path,"rice.csv") )
Group = read.csv(file.path(path, "group_rice.csv") ) [,2]

# default parameters
RCutoff = .85
MCutheight = .1
PowerUpper = 20
minModuleSize = 20

IgnoreCols = 4 	

enableWGCNAThreads()

# Only abundance data - samples as rows and proteins as columns
datExpr = as.data.frame(t(allDat[,-c(1:IgnoreCols)]) )

rownames(datExpr) = colnames(allDat)[-c(1:IgnoreCols)]
colnames(datExpr) = allDat[,1]

# Log transformation 
datExpr = log(datExpr) 	

#########################
# Build WGCNA network
#########################
		
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=PowerUpper, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = RCutoff, verbose = 5)

# Plot the results:
png("ScaleFreeTopology.png", 3000,2000,res=300)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	 xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	 main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	 labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
h = 0.85
abline(h=h,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
	 xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	 main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")	 

dev.off()

# Get the soft power
softPower = sft$powerEstimate;

# Build the adjacency table - use "signed" for proteomics data
adjacency = adjacency(datExpr, power = softPower, type="signed");


# Turn adjacency into topological overlap distance
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM


# Clustering using TOM-based dissimilarity

proTree = hclust(as.dist(dissTOM), method = "average");


# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = proTree, distM = dissTOM, 
				deepSplit = 2, pamRespectsDendro = FALSE,
				minClusterSize = minModuleSize);

print("Dynamic tree cut results:")
print(table(dynamicMods))


# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(proTree, dynamicColors, "Dynamic Tree Cut",
					dendroLabels = FALSE, hang = 0.03,
					addGuide = TRUE, guideHang = 0.05,
					main = "Protein dendrogram and module colors")

# Merge clusters

mergedClust = mergeCloseModules(datExpr, dynamicColors, cutHeight = MCutheight, verbose = 3)

mergedColors = mergedClust$colors;

mergedMEs = mergedClust$newMEs;
		 
				
	#############################################
	# Calculate and plot module eigenproteins
	# - Dendrogram
	# - Heatmap
	# - Boxplot
	# - KME
	#############################################


# Rename to moduleColors
moduleColors = mergedColors


print("Modules after merging:")
print(table(moduleColors))

# Plot dendrogram
	
png("DendroColorMergedClust.png", 2000, 2000, res=300)
plotDendroAndColors(proTree, cbind(dynamicColors, mergedColors),
					c("Dynamic Tree Cut", "Merged dynamic"),
					dendroLabels = FALSE, hang = 0.03,
					addGuide = TRUE, guideHang = 0.05)
dev.off()
	 
# Get the module eigenproteins
MEs = mergedMEs

rownames(MEs) = rownames(datExpr)


# reorder MEs by color names of modules
MEs = MEs[,order(colnames(MEs))]

# Plot module profiles with eigenproteins overlaid
WGCNAClusterID = moduleColors

plotClusterProfileWGCNA(t(datExpr), WGCNAClusterID, Group,  MEs= MEs,
		ylab="Average log ratio", file="WGCNAClusterPattenME.png",
		cex.main=1.8, cex.lab=1.7, cex.axis=1.5)

		
plotClusterProfileWGCNA(t(datExpr), WGCNAClusterID, Group,  
		ylab="Average log ratio", file="WGCNAClusterPattenAve.png",
		cex.main=1.8, cex.lab=1.7, cex.axis=1.5)


# dendrogram and heatmap for eigenproteins
png("Dendrogram eigenproteins.png", 500,500,res=100)					
plotEigengeneNetworks(MEs, "Eigenprotein Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,2,0),
	plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))	
dev.off()


png("Heatmap eigenproteins.png", 550,500,res=100)
heatmap3(t(MEs), #distfun = function(x) dist(x, method="euclidean"), 
		ColSideColors=rainbow(nlevels(as.factor(Group)))[as.factor(Group)],
		method = "average", 
		main="Module eigenproteins")
legend("topleft", fill=rainbow(nlevels(as.factor(Group)))[1:nlevels(as.factor(Group))],
			legend=levels(as.factor(Group)), cex=.6, xpd=TRUE, inset=-.1 )
dev.off()


# Network heatmap

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

png("Network heatmap.png", 2000, 2000, res=300)
TOMplot(plotTOM, proTree, moduleColors, main = "Network heatmap plot, all proteins")
dev.off()

	
# Export networks for visulisation with VisANT - can take a while!
if(FALSE) {	
	
	for(module in moduleColors) {
		# Select module probes
		probes = names(datExpr)
		inModule = (moduleColors==module);
		modProbes = probes[inModule];
		
		# Select the corresponding Topological Overlap
		
		for(NmodTOM in 1:2) {			
			
			modTOM = TOM[inModule, inModule];
			visFile = paste("VisANTInput-TOM", module, ".txt", sep="")
			
			if(NmodTOM == 2) {
				modTOM = adjacency[inModule, inModule]
				visFile = paste("VisANTInput-ADJ", module, ".txt", sep="")
			}	

			dimnames(modTOM) = list(modProbes, modProbes)
			# Export the network into an edge list file VisANT can read
			vis = exportNetworkToVisANT(modTOM,
				  file = visFile,
				  weighted = TRUE,
				  threshold = 0
			  )			
		}
	}
}

# Get KME - module membership - correlation between proteins and eigenproteins
	
kmes = signedKME(datExpr, MEs)


# separate results by modules, order by kME, hub proteins on top

dat.res = data.frame(allDat, moduleColors , kmes)

list.cluster.dat = lapply(levels(as.factor(moduleColors)), 
	function(x) {dtemp = dat.res[dat.res$moduleColors == x,];
			dtemp[order(dtemp[,paste0('kME',x)==colnames(dtemp)], decreasing=TRUE),
				-setdiff(grep("^kME", colnames(dtemp)), which(paste0('kME',x)==colnames(dtemp)))]} )
			
names(list.cluster.dat) = 	levels(as.factor(moduleColors))



# Boxplot for eigenproteins

ag.temp = aggregate(MEs, by=list(Group=Group), FUN=mean)
ag.eigengenes = t(ag.temp[,-1])
colnames(ag.eigengenes) = ag.temp[,1]

fScale = max(8,nlevels(as.factor(moduleColors)))/8

png("Boxplot eigenproteins.png", 2000, 3000*fScale, res=300)

par(mar= c(7, 4, 2, 1) + 0.1)
layout(matrix(1:(ceiling(nlevels(as.factor(moduleColors))/2)*2), ncol=2, byrow=TRUE))
cols = levels(as.factor(moduleColors))
for(ii in 1:ncol(MEs))	
	boxplot(MEs[,ii] ~ Group, las=2, col=cols[ii], ylab = "log ratio",
	main=paste(colnames(MEs)[ii], table(moduleColors)[ii] ), cex.main=1.7, cex.lab=1.7, cex.axis=1.5 )

dev.off()

# Boxplot for top 6 hub proteins

for(ii in 1:length(list.cluster.dat)) {

png(paste0("Boxplot hub proteins - ", names(list.cluster.dat)[ii], ".png"), 2000, 2500, res=300)
par(oma= c(5, 2, 2, 1) + 0.1)
layout(matrix(1:6, ncol=2))
for(jj in 1:6) boxplot(t(log(list.cluster.dat[[ii]][jj,5:22])) ~ Group, 
	main=paste(list.cluster.dat[[ii]][jj,2],"\nkME=", round(list.cluster.dat[[ii]][jj,ncol(list.cluster.dat[[ii]])],2)), 
	col=rainbow(nlevels(as.factor(Group))), ylab="Log ratio", cex.main=1.5, las=2,cex.lab=1.2, )


dev.off()

}

				
# Output results
		
wb = createWorkbook()

addWorksheet(wb, "AllData")

writeData(wb, "AllData", dat.res)

# write modules only tabs

for(ii in 1:length(list.cluster.dat)) {
	addWorksheet(wb, names(list.cluster.dat)[ii])
	writeData(wb, names(list.cluster.dat)[ii], list.cluster.dat[[ii]])

}

saveWorkbook(wb, "ResultsWGCNA.xlsx", overwrite=TRUE)

# output the eigenprotein	
write.csv(MEs,"Module eigenprotein.csv")

