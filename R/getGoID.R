getGoID <- function(v) {

xx <- as.list(GOTERM)
tp <- sapply(xx, FUN=Term)
tp[match(tolower(v), tolower(tp), nomatch=0)]

}

