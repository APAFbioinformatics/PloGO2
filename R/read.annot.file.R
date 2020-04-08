read.annot.file <- function(fname, format=c("compact","long") ) {

    format <- match.arg(format)
###############################
# accept space or ; separators
###############################

tbl <- read.table(fname, sep = "\n")
IDS <- apply(as.matrix(tbl), 1, FUN = function(v) {strsplit(as.vector(v), ";|[[:space:]]+")[[1]][1] })

if (format == "long") {
       # consider the second value only, discard the rest
        V1 <- apply(as.matrix(tbl), 1, FUN = function(v) {strsplit(as.vector(v), "\t")[[1]][2] })
        tbl <- as.matrix(data.frame(IDS,V1))
        annot.dat <- aggregate(tbl[, 2], by = list(IDS = tbl[,1]), FUN = function(x) {paste(x, collapse = " ") })
        colnames(annot.dat) <- c("IDS", "V1")
} else {
        # consider all remaining values
        V1 <- apply(as.matrix(tbl), 1, FUN = function(v) { paste(strsplit(as.vector(v), ";|[[:space:]]+")[[1]][-1], collapse = " ") })
        annot.dat <- data.frame(IDS, V1)
}
    
annot.dat



}
