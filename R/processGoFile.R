processGoFile <- function(fname, GOIDlist, datafile=NULL, datafile.ignore.cols=1, 
	format="compact", aggregateFun="sum") {


termListIdx <- list()
annot.dat <- read.annot.file(fname, format=format)
abundance  <- NULL
N <- nrow(annot.dat)

# for each term, extract an index: which proteins belong to it
for (term in GOIDlist) {
  g1 <- GOGraphWrapper(term)
  term.idx <- apply(as.matrix(annot.dat[,2]),1, FUN=function(v){inGraph(v,g1)})
  termListIdx[[term]] <- term.idx
}

# extract the actual proteins
full.ID.list <- lapply(termListIdx, FUN=function(v){as.vector(annot.dat[v,1])})


counts <- sapply(termListIdx, FUN=sum)

names <- sapply(names(termListIdx), FUN=function(v){Term(GOTERM[[v]])})
counts <- data.frame(counts, names)
names(counts)[1] <- fname



if (!is.null(datafile)) {

nsaf <- read.csv(datafile)

if (aggregateFun == "prod") {

abundance <- sapply(full.ID.list, FUN=function(v){apply(nsaf[match(unlist(v), nsaf[,1], nomatch=0) ,-c(1:datafile.ignore.cols)],2,FUN=prod)})

} else {

abundance <- sapply(full.ID.list, FUN=function(v){apply(nsaf[match(unlist(v), nsaf[,1], nomatch=0) ,-c(1:datafile.ignore.cols)],2,FUN=sum)})

}

} 


list(counts=counts, ID.list=full.ID.list, datafile=datafile, abundance=abundance, N=N, fname=gsub(".*\\/", "", fname)) 


}


# res <- processGoFile(fname,GOIDlist)
