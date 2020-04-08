
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
    # 
    list.barplots <- list()
    
     if(Plot) {
    for (gg in unique(GG)[seq_len(min(length(unique(GG)),printLimit)) ] ) {
      # #JW
    				
	list.barplots[[gg]] <- barchart(Value ~ GG | FF, groups = SS, , main = gg,
            scales = list(x = list(rot = 45)), auto.key = list(points = FALSE,
                rectangles = TRUE, space = "top"), subset = (GG ==
                 gg))

    }

    }

    # Only plot counts >= CountCutOff
    list.levelplots <- list()
    
    for (i in seq_along(fnames)) {
        
        nzero = sum(!(counts.list[[i]]==0))
        idx2plot = counts.list[[i]] >= CountCutOff
        
        if(sum(idx2plot) > 2) {            

          # png(paste("File", i, ".png", sep = ""), 2000, 4000, res = 200)
       
       list.levelplots[[fnames[i]]] <- plotMat(abundance.list[[i]][,idx2plot], counts.list[[i]][idx2plot], 
			main = paste("Abundance levelplot",
            fnames[i]))
        # dev.off()
        
        }
    }

	ag.mat <- NULL
	
    if (!is.null(Group)) {

	ag.list <- aggregateAbundance(res.list, Group=Group)
	#
	# should review this part; DP Nov 2019
	# remove 
	x <- extractPairs(res.list, ag.list, fnames[1:2], Group[1])
	
	# png("CountAndAbundance.png", 2000, 2000, res=300)
	# countAndAbundance(x, legend.text=fnames[1:2])
	# dev.off()

	ag.mat <- ag.list[[1]]
	
	for (jj in 2:length(ag.list)) ag.mat <- rbind(ag.mat, ag.list[[jj]]);
	
	Files <- as.factor(rep(fnames, each=nlevels(as.factor(Group))))
	rownames(ag.mat) <- paste(rownames(ag.mat), Files)

	#write.csv(ag.mat, file="AggregateAbundanceMat.csv")

	}

	 # abundance
   list(abundance=abundance, ag.mat=ag.mat, list.levelplots=list.levelplots, list.barplots=list.barplots)
}
