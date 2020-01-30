inferGroup <- function(v) {

v <- gsub("^([\\. ])*", "",  v)

v <- gsub("\\.csv$", "", v)

v <- gsub("[[:digit:]]*$", "", v)

v <- tolower(gsub("\\.*$", "", v))

as.factor(v)

}
