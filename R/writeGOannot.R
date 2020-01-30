writeGOannot <- function(res, fname="AnnotOut.txt", datafile = NULL, format="list") {


# two formats: list, or adjacency matrix 


ID.list <- res$ID.list
counts <- as.matrix(res$counts)

if (!is.null(datafile)) additionalData <- read.csv(datafile);

if (format == "matrix") {

# output to adjacency matrix
   IDS <- unique(unlist(ID.list))
        adj.mat <- sapply(ID.list, FUN = function(v) {
            w <- rep(0, length(IDS))
            w[match(v, IDS)] <- 1
            w
        })
        adj.mat <- data.frame(IDS, adj.mat)
	descriptions <- c(NA,counts[,2])
        if (!is.null(datafile)) {
            adj.mat <- merge(adj.mat, additionalData, by.x = 1, 
                by.y = 1, all.x = TRUE)
	    dataDesc <-  colnames(additionalData)[-1]
	    descriptions <- c(descriptions,dataDesc )
        }

        adj.mat <- data.frame( descriptions , t(adj.mat))
        write.csv(t(adj.mat), file = fname, row.names = FALSE)

} else {

if (!is.null(datafile)) {

ID.list <- lapply(ID.list, FUN=function(v){merge(v, additionalData, by.x=1, by.y=1, all.x=TRUE)})

}


cat("Annotation summary by category \n", file=fname)

for (i in 1:length(ID.list)) {
	cat(rownames(counts)[i], counts[i,2], "\n", 
		"Number:", counts[i,1], "\n", sep="\t", 
		file=fname, append=TRUE)
	write.table(format(ID.list[[i]], digits=4), sep="\t", file=fname, append=TRUE)

}


}  # end output to text file

}


# writeGOannot(res, datafile=datafile)
