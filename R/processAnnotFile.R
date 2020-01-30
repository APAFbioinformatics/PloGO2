processAnnotFile <- function (fname, GOIDlist, datafile = NULL, datafile.ignore.cols = 1, 
    format = "compact", term.names=NULL) 
{
    termListIdx <- list()
    annot.dat <- read.annot.file(fname, format = format)
    abundance <- NULL
    N <- nrow(annot.dat)
    for (term in GOIDlist) {
        term.idx <- apply(as.matrix(annot.dat[, 2]), 1, FUN = function(v) {
            v == term
        })
        termListIdx[[term]] <- term.idx
    }
    full.ID.list <- lapply(termListIdx, FUN = function(v) {
        as.vector(annot.dat[v, 1])
    })
    counts <- sapply(termListIdx, FUN = sum)




   if (is.null(term.names)) {
	names <- names(termListIdx)

	}  else {
    # validate term names?
    names <- sapply(names(termListIdx), FUN = function(v) {
        as.vector(term.names[match(v, term.names[,1]),2])
    })
    }

    counts <- data.frame(counts, names)
    names(counts)[1] <- fname
    if (!is.null(datafile)) {
        nsaf <- read.csv(datafile)
        abundance <- sapply(full.ID.list, FUN = function(v) {
            apply(nsaf[match(unlist(v), nsaf[, 1], nomatch = 0), 
                -c(1:datafile.ignore.cols)], 2, FUN = sum)
        })
    }
    list(counts = counts, ID.list = full.ID.list, datafile = datafile, 
        abundance = abundance, N = N, fname = gsub(".*\\/", "", 
            fname))
}

