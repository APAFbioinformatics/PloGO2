annotationPlot <- function (res.list, percentages = FALSE, plot = TRUE, trimzero = FALSE,
		type=c("GO", "pathway") )
{
	type = match.arg(type)
   counts <- sapply(res.list, FUN = function(v) {
        v$counts[, 1]
    })
    N <- sapply(res.list, FUN = function(v) {
        v$N
    })
    fnames <- sapply(res.list, FUN = function(v) {
        v$fname
    })
    colnames(counts) <- fnames
    Percentages <- 100 * sweep(counts, 2, N, FUN = "/")
    GOIDlist <- rownames(res.list[[1]]$counts)
	goTerms = GOIDlist
	if(type == "GO")
    goTerms <- sapply(GOIDlist, FUN = function(v) {
        Term(GOTERM[[v]])
    })
    counts.print <- counts
    rownames(counts.print) <- goTerms
    Percentages.print <- Percentages
    rownames(Percentages.print) <- goTerms

	if (trimzero) {
		z.idx <- (apply(counts.print, 1, FUN=sum) == 0)
		counts.print <- counts.print[!z.idx, ]
		Percentages.print <- Percentages.print[!z.idx, ]
	}

    #write.csv(counts.print, file = "Counts.csv")
    #write.csv(Percentages.print, file = "Percentages.csv")
    rownames(counts) <- GOIDlist
    rownames(Percentages) <- GOIDlist
    AnnotNumbers <- as.vector(counts)
    AnnotPercentages <- as.vector(Percentages)
    GOID <- rep(GOIDlist, ncol(counts))
    GOTerm <- rep(goTerms, ncol(counts))
    Files <- rep(fnames, each = nrow(counts))
    df <- data.frame(AnnotNumbers, AnnotPercentages, GOID, GOTerm, Files)

    list(counts = counts.print, percentages = Percentages.print)
}
