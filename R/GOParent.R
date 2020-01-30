
GOParent <- function (node) 
{

parent <- NULL

if (!is.null(GOTERM[[node]])) {
    res <- get(paste("GO", "PARENTS", sep = Ontology(GOTERM[[node]])))[[node]]
    # depending on GOstats version labels might differ as below
    # for "is a" relationship in GO graph
    parent <- res[names(res) %in% c("is.a", "is_a", "isa")][1]
}

parent

}

