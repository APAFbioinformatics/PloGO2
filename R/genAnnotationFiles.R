
# Generate a folder of Wengo format annotation files from an excel input and a DB file

# Example: genAnnotationFiles("allStrains.xlsx")

library(openxlsx)

genAnnotationFiles <- function(fExcelName,  colName="Uniprot", 
				DB.name = "pathwayDB.csv", folder="PWFiles"
				)
{
	
	sheets = names(loadWorkbook(fExcelName))
	
	data.list = lapply(sheets, FUN=function(s){ readWorkbook(fExcelName, s) } )
	id.list <- lapply(data.list, FUN=function(d){unique(d[,colName])})
	# generate GO files 
	if ( !( folder %in% list.files()) ) dir.create(folder)

	# infer identifier type for Uniprot or Ensembl only!
	IDS <- unique(unlist(id.list))

	# Read in DB file
	dat_db = read.csv(DB.name, stringsAsFactors=FALSE)
	
	
	for (ii in 1:length(id.list) ) {
		vvv <- unique(id.list[[ii]])
		
		ANNOTID <- rep("", length(vvv))
		ANNOTID[!is.na(match(vvv, dat_db[,1]))] <- dat_db[match(vvv, dat_db[,1], nomatch = 0),2]
		res <- data.frame(vvv, ANNOTID)
		names(res) <- c("Identifier", "ID")
		write.table(res, file = paste0(folder,"/",sheets[ii],".txt"), sep = "\t", quote = FALSE, 
			row.names = FALSE, col.names = FALSE)
		
	}


}
				