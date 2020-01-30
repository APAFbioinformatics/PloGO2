GOGraphWrapper <- function(term) {

	GOGraph(term, get(paste("GO", "CHILDREN", sep = Ontology(GOTERM[[term]])))) 

}
