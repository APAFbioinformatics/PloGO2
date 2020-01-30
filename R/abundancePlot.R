
abundancePlot <- function (res.list, log = FALSE, printLimit = 16, Group=NULL, Plot=FALSE, 
	CountCutOff=3, ...)
{
    abundance.list <- lapply(res.list, FUN = function(v) {
        v$abundance
    })
    counts.list <- lapply(res.list, FUN = function(v) {
        v$counts[, 1]
    })
    abundance <- sapply(res.list, FUN = function(v) {
        v$abundance
    })
    fnames <- sapply(res.list, FUN = function(v) {
        v$fname
    })
    colnames(abundance) <- fnames
    
    IDlist <- colnames(res.list[[1]]$abundance)
  

# DP 2012-05-24
# modify to allow for random lists, not just GO

    if (length(grep("^GO:", IDlist)) == length(IDlist)) {
      #GOIDlist = gsub(";", "", GOIDlist)
      
	Terms <- sapply(IDlist, FUN = function(v) { Term(GOTERM[[v]]) } )
	#jw 11/11/19
	IDlist <- Terms
    } else {
	Terms <- IDlist
    }


    nFiles <- rownames(res.list[[1]]$abundance)
    Source <- rep(nFiles, length(IDlist))
    Term <- rep(Terms, each = length(nFiles))
    rownames(abundance) <- paste(Source, Term)
    Value <- as.vector(abundance)
    if (log) Value <- log(Value + eps);
    GG <- rep(Term, ncol(abundance))
    SS <- rep(Source, ncol(abundance))
    FF <- rep(gsub(".*\\/", "", colnames(abundance)), each = nrow(abundance))   
    
    if(Plot) {
    for (gg in unique(GG)[1:min(length(unique(GG)),printLimit) ] ) {
      #JW
      fImage = gsub(" ", "", paste(substring(gg, 1, 30), ".png", sep = ""))
      fImage = gsub("/", " ", fImage)
      
        # png(gsub(" ", "", paste("Categ", substring(gg, 1, 30), ".png", sep = "")),
      png(fImage, 2000, 2000, res = 200)
        print(barchart(Value ~ GG | FF, groups = SS, , main = gg,
            scales = list(x = list(rot = 45)), auto.key = list(points = FALSE,
                rectangles = TRUE, space = "top"), subset = (GG ==
                gg)))
        dev.off()
    }
    
    }

    # Only plot counts >= CountCutOff
    for (i in 1:length(fnames)) {
        
        nzero = sum(!(counts.list[[i]]==0))
        idx2plot = counts.list[[i]] >= CountCutOff
        
        if(sum(idx2plot) > 2) {            

          png(paste("File", i, ".png", sep = ""), 2000, 4000, res = 200)
       
        plotMat(abundance.list[[i]][,idx2plot], counts.list[[i]][idx2plot], 
			main = paste("Abundance levelplot",
            fnames[i]), log=log, ...)
        dev.off()
        
        }
    }

	
    if (!is.null(Group)) {

	ag.list <- aggregateAbundance(res.list, Group=Group)
	#
	# should review this part; DP Nov 2019
	# remove 
	x <- extractPairs(res.list, ag.list, fnames[1:2], Group[1])
	
	png("CountAndAbundance.png", 2000, 2000, res=300)
	countAndAbundance(x, legend.text=fnames[1:2])
	dev.off()

	ag.mat <- ag.list[[1]]
	
	for (jj in 2:length(ag.list)) ag.mat <- rbind(ag.mat, ag.list[[jj]]);
	
	Files <- as.factor(rep(fnames, each=nlevels(as.factor(Group))))
	rownames(ag.mat) <- paste(rownames(ag.mat), Files)

	write.csv(ag.mat, file="AggregateAbundanceMat.csv")

	#
	# Abundance barplot DP JW Nov 2019
	#
	# a new function now plotAbundanceBar
	if(Plot) {
	ag.mat[ag.mat==0] = NA		 

	Objects = IDlist
	counts <- sapply(res.list, FUN = function(v) v$counts[, 1])
	
	Max = apply(counts, 1, FUN=max)

	keep.idx = (Max > 4) #( Max > 3 )
	obj = Objects[keep.idx]

	
	ag.mat = ag.mat[, (tolower(colnames(ag.mat)) %in% tolower(obj)), drop=FALSE]
	N = ncol(ag.mat)
	
	png("Abundance barplot.png", 2500, 2000, res=300)
	par(mar=c(4,10,4,14))
	
	barplot(as.matrix(t(ag.mat)), las=2, col=rainbow(N), horiz=TRUE, cex.names=.6)
	legend("right", fill=rainbow(N), legend=colnames(ag.mat), cex=0.6,xpd=TRUE, inset=-.6)

	dev.off()
  }
	}




	 abundance
}
