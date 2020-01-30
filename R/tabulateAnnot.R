

tabulateAnnot <- function(res.list) {

ID.list <- res.list$ID.list
counts <- as.matrix(res.list$counts)
if(length(grep("^GO:", names(ID.list))) == length(ID.list)) {
  Names <- unlist(lapply(names(ID.list), FUN=function(x){Term(GOTERM[[x]])}))
} else {
  Names <- rownames(counts)[match(names(ID.list), counts[,3])]
}

names(ID.list) <- paste(Names, names(ID.list))
all.ID <- unique(unlist(ID.list))
adj.mat <- sapply(ID.list, FUN = function(v) {
            w <- rep(0, length(all.ID))
            w[match(v, all.ID)] <- 1
            w
})

data.frame(all.ID, adj.mat)

}
