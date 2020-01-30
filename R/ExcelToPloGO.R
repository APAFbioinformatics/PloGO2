

ExcelToPloGO <- function(fname,  colName="Uniprot", 
				termFile=NA,
				compareWithReference="none",
				data.file.name = "none"
				)

{

# if no files uploaded, default categories
if (is.na(termFile)) {
 path <- system.file("files", package = "PloGO2")
 termFile = paste(path, "GODefault.txt", sep="/")
}

wb <- (fname)
# sheets <- getSheets(wb)

sheets <- names(openxlsx::loadWorkbook(wb))

data.list = lapply(sheets, FUN=function(s){ openxlsx::readWorkbook(wb, s) } )
id.list <- lapply(data.list, FUN=function(d){unique(d[,colName])})
# generate GO files 
if ( !( "./GOFiles" %in% list.files()) ) dir.create("GOFiles")

# infer identifier type for Uniprot or Ensembl only!
IDS <- unique(unlist(id.list))
database <- "uniprot"
if ( length(grep("^ensp", tolower(IDS))/max(1,length(IDS)))  > .5 ) database <- "ensembl";

for (ii in 1:length(id.list) ) {
vvv <- unique(id.list[[ii]])
genWegoFile(vvv, database=database,fname=paste("GOFiles/", sheets[[ii]], ".txt", sep=""))
}

plogo.res <- PloGO( filesPath="GOFiles", termFile=termFile, 
				reference=compareWithReference, data.file.name=data.file.name )

# merge counts, percentages and Pvalues
# add GO category


Counts = plogo.res$Counts
colnames(Counts) = paste("Counts", colnames(Counts))
Percentages = plogo.res$Percentages
colnames(Percentages) = paste("Percentages", colnames(Percentages))
results = data.frame(Counts, Percentages)
if ( !is.null(plogo.res$FisherPval) ){
FisherPval = plogo.res$FisherPval[rownames(plogo.res$FisherPval) %in% rownames(plogo.res$Counts),]
colnames(FisherPval) = paste("FisherPval", colnames(FisherPval))
results = data.frame(results, FisherPval)
}

xx <- as.list(GOTERM)
tp <- sapply(xx, FUN = Term)
Ont <- sapply(xx, FUN = Ontology)
match.idx = match(tolower(rownames(results)), tolower(tp), nomatch = 0)
go.annot = data.frame(tp[match.idx], Ont[match.idx])
colnames(go.annot) = c("GO.Term", "Ontology")

results <- data.frame(go.annot, results)


annotated.list <- list()

# make sure the folder order of files from PloGO is the
# same as the sheet order from the spreadsheet!
# sheet.ord <- match(paste(sheets, ".txt", sep=""), names(plogo.list[[1]]$res.list))

sheet.ord <- match( names(plogo.res$res.list), paste(sheets, ".txt", sep=""))
data.list <- data.list[sheet.ord]

for ( iii in 1:length(data.list) ) {

dat.res <- data.list[[iii]]
dat.res$Order <- 1:nrow(dat.res)
adj.mat <- tabulateAnnot( plogo.res$res.list[[iii]] )
dat.res <- merge( dat.res, adj.mat, by.x=colName, by.y=1, all.x=TRUE, sort=FALSE)
annotated.list[[iii]] <- dat.res[order(dat.res$Order),]

}

# Print results to tabbed sheet
names(annotated.list) <- sheets[sheet.ord]


wb <- openxlsx::createWorkbook("Results.xlsx")


 for (jj in 1:length(annotated.list)) {
  if ("StouffersPval" %in% names(annotated.list[[jj]])) {
                pvals <- c(grep("Stouffer", names(annotated.list[[jj]])), 
                  grep("\\.pval\\.1..\\.1..$", tolower(colnames(annotated.list[[jj]]))))
                ratios <- c(grep("Geomean", names(annotated.list[[jj]])), 
                  grep("\\.X1..\\.1..$", colnames(annotated.list[[jj]])))
                ratios <- setdiff(ratios, grep("^Crit", colnames(annotated.list[[jj]])))
                tabName <- names(annotated.list)[jj]
            }
           else {
                ratios <- grep("^X1..\\.1..$", colnames(annotated.list[[jj]]))
                pvals <- grep("^pval\\.1..\\.1..$", tolower(colnames(annotated.list[[jj]])))
                tabName <- names(annotated.list)[jj]
            }
            printOpenxlsxStyle(annotated.list[[jj]], ratios, pvals, wb = wb, 
                tabName = tabName)
  }


  tabName = "GOSummary"

printSimpleOpenxlsxStyle(results, "GOSummary", wb = wb)


openxlsx::saveWorkbook(wb, file="Results.xlsx",  overwrite=TRUE)




plogo.res

}

