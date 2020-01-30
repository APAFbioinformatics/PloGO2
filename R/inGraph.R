inGraph <- function(terms, graph) {

  t <- strsplit(terms, ";|[[:space:]]+")[[1]]
  sum(nodes(graph) %in% t) > 0

}
