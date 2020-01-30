
# plot abundance by category and file
plotAbundance <- function(res.list, log=FALSE, eps = 0.001) {

    nFiles <- rownames(res.list[[1]]$abundance)
	
	IDlist <- colnames(res.list[[1]]$abundance)


    if (length(grep("^GO:", IDlist)) == length(IDlist)) {
      
	Terms <- sapply(IDlist, FUN = function(v) { Term(GOTERM[[v]]) } )
    } else {
	Terms <- IDlist
    }
	
	abundance <- sapply(res.list, FUN = function(v) {
        v$abundance
    })
	
	fnames <- sapply(res.list, FUN = function(v) {
        v$fname
    })
	
    colnames(abundance) <- fnames
	
    Source <- rep(nFiles, length(IDlist))
    Term <- rep(Terms, each = length(nFiles))
    rownames(abundance) <- paste(Source, Term)
    Value <- as.vector(abundance)
	
    if (log) Value <- log(Value + eps);
    GG <- rep(Term, ncol(abundance))
    SS <- rep(Source, ncol(abundance))
    FF <- rep(gsub(".*\\/", "", colnames(abundance)), each = nrow(abundance))    


    png("AbundanceByCategory.png", 2000, 2000, res = 300)
    print(barchart(Value ~ GG, groups = FF, scales = list(x = list(rot = 45)),
        auto.key = list(points = FALSE, rectangles = TRUE, space = "top"),
        main = "Abundance by category"))
    dev.off()

	png("AbundanceByFile.png", 2000, 2000, res = 300)
    print(barchart(Value ~ FF, groups = GG, scales = list(x = list(rot = 45)),
        auto.key = list(points = FALSE, rectangles = TRUE, space = "top"),
        main = "Abundance by file"))
    dev.off()
    
}