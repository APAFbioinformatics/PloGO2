

GOTermList <- function (ontology = "BP", level = 2, node = NULL) 
{
    root.terms <- list(MF = "GO:0003674", CC = "GO:0005575", 
        BP = "GO:0008150")
    if (!is.null(node)) {
        ontology <- Ontology(GOTERM[[node]])
        res <- c(node, Children(node))
    }
    else {
        root <- root.terms[[ontology]]
        if (level == 1) 
            res <- root
        if (level == 2) 
            res <- Children(root)
        if (level == 3) {
            lev2 <- Children(root)
	    intermed <- lapply(lev2, FUN = Children)
	    tp <- unlist(intermed)
	    # extract unique; 
	    uniq <- names(table(tp))
	    res <- tp[match(uniq, tp)]
        }
        if (level == 4) {
	    lev2 <- Children(root)
	    intermed <- lapply(lev2, FUN = Children)
	    # to ensure the right names
	    tp <- unlist(intermed)
	    # extract unique; retain one parent
	    uniq <- names(table(tp))
	    lev3 <- tp[match(uniq, tp)]
	    intermed <- lapply(lev3, FUN = Children)
	    tp <- unlist(intermed)
	    uniq <- names(table(tp))
	    res <- tp[match(uniq, tp)]

        }
    }
    res
}


