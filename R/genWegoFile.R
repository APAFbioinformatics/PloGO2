
genWegoFile <- function (IDList, fname = "Wego.txt", database="uniprot") 
{


# remove duplicates
IDList = IDList[!duplicated(IDList)]


CleanUniprot <- function(vec){
Uniprot = gsub("(sp|tr)\\|(.*)\\|(.*)", "\\2", vec)
Uniprot =  gsub("^1\\/", "", Uniprot)  # clean up extended
Uniprot
}



if (database == "uniprot") {

# trembl? other IDs? Ensembl?
Uniprot = CleanUniprot(IDList)

# remove empty, NA entries
Uniprot = Uniprot[!is.na(Uniprot)]
Uniprot = Uniprot[Uniprot != ""]


# Break in batches of 250
from = seq(1,length(Uniprot), 250)
to = c(from[-1], 1+length(Uniprot)) -1

unip.list = list()
for (jj in c(1:length(from))) {
   # unip.list[[jj]] = try( getUniprotBatch(values=Uniprot[from[jj]:to[jj]], attributes=c("id", "protein+names", "go")) )
   unip.list[[jj]] = try( getUniprotBatch(values=Uniprot[from[jj]:to[jj]], attributes=c("id", "go-id")) )
} 

for (jj in 1:length(from)) {
if (jj == 1) { annotated = unip.list[[1]] } else { annotated = rbind(annotated, unip.list[[jj]]) } ;
}

bm.res = annotated
	
		
} else {
    mart <- useMart("ensembl")
    mart<-useDataset("hsapiens_gene_ensembl",mart)


bm <- getBM(filters = "ensembl_peptide_id", values = IDList, 
            attributes = c("ensembl_peptide_id", "go_id"), 
            mart)
bm.res <- aggregate(bm[, 2], by = list(Acc = bm[, 1]), FUN = function(v){paste(v, collapse = " ")})
   
}

    GOID <- rep("", length(IDList))
    GOID[match(bm.res[, 1], IDList, nomatch = 0)] <- bm.res[, 
        2]
    res <- data.frame(IDList, GOID)
    names(res) <- c("Identifier", "ID")
    write.table(res, file = fname, sep = "\t", quote = FALSE, 
        row.names = FALSE, col.names = FALSE)
		
		
}



