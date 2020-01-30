
aggregateAbundance <- function(res.list, Group) {



ag.list <- list()

if ( (length(res.list) > 0) & (!is.null(res.list[[1]]$abundance)) ) {

if (nrow(res.list[[1]]$abundance) != length(Group) ) Error("The Group length does not match the number of data columns of the annotation matrix");
 


GOIDlist <- colnames(res.list[[1]]$abundance)

counts <- sapply(res.list, FUN = function(v) {
        v$counts[, 1]
    })

fnames <- sapply(res.list, FUN = function(v) {
        v$fname
    })


# DP 2012-05-24
# modify to allow for random lists, not just GO

if (length(grep("^GO:", GOIDlist)) == length(GOIDlist)) {
	goTerms <- sapply(GOIDlist, FUN = function(v) {
        Term(GOTERM[[v]])
    })
} else {
	goTerms <- GOIDlist
}



rownames(counts) <- goTerms
colnames(counts) <- fnames

abundance.list <- lapply(res.list, FUN = function(v) {
         v$abundance
     })

names(abundance.list) <- fnames

ag.list <-  lapply(abundance.list, FUN=function(x){
                mat <- aggregate(x, by=list(Group=Group), FUN=mean)
                rownames(mat) <- mat[,1]
                mat <- mat[,-1]
                colnames(mat) <- goTerms
                mat
                })

names(ag.list) <- fnames

}


ag.list



}
