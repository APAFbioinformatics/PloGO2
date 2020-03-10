

ExcelToPloPathway <- function(fname,  colName="Uniprot", 
				compareWithReference="none", DB.name="pathwayDB.csv",
				data.file.name = "none")
{

	wb <- (fname)

	sheets <- names(openxlsx::loadWorkbook(wb))

	data.list = lapply(sheets, FUN=function(s){ openxlsx::readWorkbook(wb, s) } )
	id.list <- lapply(data.list, FUN=function(d){unique(d[,colName])})

	# get the annotation files	
	genAnnotationFiles(fname,  
					colName=colName, 
					DB.name = DB.name, folder="PWFiles")

	# Run PloPathway
	res = PloPathway(zipFile = "none",
			reference = compareWithReference,
			data.file.name = data.file.name, 
			datafile.ignore.cols = 1, 
			filesPath="PWFiles")

	if(!data.file.name=="none")		
	# Plot abundance barplot
	plotAbundanceBar(res$aggAbundance[-c(1:4),], res$Counts[,-1], min.count=4)


	# Print results	in temp file
	printSummary(res, file="PloGO2Results.xlsx")		

	sheets.plogo = names(loadWorkbook("PloGO2Results.xlsx"))
	
	if(!length(sheets) == (length(sheets.plogo)-1) ) warnings("the number of sheets in the input is different from PloGO2Results.xlsx")

	# merge plopathway results and input file
	list.mg = lapply(1:length(sheets), function(ii) 
		merge(readWorkbook("PloGO2Results.xlsx", ii), data.list[[ii]], by.x="Protein", by.y=colName, 
		all=TRUE, sort=FALSE))

	wb <- createWorkbook()

	for(ii in 1:length(list.mg)) {
		addWorksheet(wb, sheets[ii])
		writeData(wb, sheets[ii], list.mg[[ii]])
		
	}

	addWorksheet(wb, sheets.plogo[length(sheets.plogo)])
	writeData(wb, sheets.plogo[length(sheets.plogo)],
		readWorkbook("PloGO2Results.xlsx", length(sheets.plogo)) )

	saveWorkbook(wb, "PloGO2Results.xlsx", overwrite=TRUE)

	res
}
				