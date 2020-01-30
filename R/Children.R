Children <- function (node) 
{

Child <- NULL

if (!is.null(GOTERM[[node]])) {
	Child <- get(paste("GO", "CHILDREN", sep = Ontology(GOTERM[[node]])))[[node]]
}

names(Child) <- NULL

Child
}
