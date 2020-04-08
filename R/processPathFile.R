

processPathFile <- function(fname, AnnotIDlist, datafile=NULL, datafile.ignore.cols=1, 
	format=c("compact","long"), aggregateFun="sum") {
format <- match.arg(format)
termListIdx <- list()
annot.dat <- read.annot.file(fname, format=format)
abundance  <- NULL
N <- nrow(annot.dat)

# for each term, extract an index: which proteins belong to it
for (term in AnnotIDlist) {
    term.idx = grepl(term, annot.dat[,2])
  
  termListIdx[[term]] <- term.idx
}

# extract the actual proteins
full.ID.list <- lapply(termListIdx, FUN=function(v){as.vector(annot.dat[v,1])})


counts <- sapply(termListIdx, FUN=sum)

names <- names(termListIdx)

# get the numeric pathway IDs
pathIDs <- as.numeric(gsub("[a-z]", "", names))
  
names = keggPathway[match(pathIDs, keggPathway[,1]),2]
  
counts <- data.frame(counts, names)
names(counts)[1] <- fname


if (!is.null(datafile)) {

datAbun <- read.csv(datafile)

if (aggregateFun == "prod") {


abundance <- sapply(full.ID.list, FUN=function(v){
	matAbun = datAbun[match(unlist(v), datAbun[,1], nomatch=0) ,-seq_len(datafile.ignore.cols)];
	for(ii in seq_len(ncol(matAbun))) matAbun[,ii] = as.numeric(matAbun[,ii]);
	apply(matAbun,2,FUN=prod)})

} else {

abundance <- sapply(full.ID.list, FUN=function(v){
	matAbun = datAbun[match(unlist(v), datAbun[,1], nomatch=0) ,-seq_len(datafile.ignore.cols)];
	for(ii in seq_len(ncol(matAbun))) matAbun[,ii] = as.numeric(matAbun[,ii]);
	apply(matAbun,2,FUN=sum)})

}
colnames(abundance) = names
} 


counts$pathwayID = rownames(counts)
rownames(counts) = counts[,2]

list(counts=counts, ID.list=full.ID.list, datafile=datafile, abundance=abundance, N=N, fname=gsub(".*\\/", "", fname)) 


}


