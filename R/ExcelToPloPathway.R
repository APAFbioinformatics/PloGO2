

ExcelToPloPathway <- function(fname,  colName="Uniprot", 
				compareWithReference="none", DB.name="pathwayDB.csv",
				data.file.name = "none",outFolder = "PloGO2Output")
{

	wb <- (fname)

	sheets <- names(openxlsx::loadWorkbook(wb))

	data.list = lapply(sheets, FUN=function(s){ openxlsx::readWorkbook(wb, s) } )
	id.list <- lapply(data.list, FUN=function(d){unique(d[,colName])})

	
	# make output folder
	if(!(outFolder %in% list.files())) dir.create(outFolder)
	
	# get the annotation files	
	filespath <- genAnnotationFiles(fname,  
					colName=colName, 
					DB.name = DB.name, folder="PWFiles", outFolder)

	# Run PloPathway
	res = PloPathway(zipFile = "none",
			reference = compareWithReference,
			data.file.name = data.file.name, 
			datafile.ignore.cols = 1, 
			filesPath=filespath)


	# Print results	in temp file
	printSummary(res, file=file.path(outFolder,"PloGO2Results.xlsx") )		

	sheets.plogo = names(loadWorkbook(file.path(outFolder,"PloGO2Results.xlsx")) )
	
	if(!length(sheets) == (length(sheets.plogo)-1) ) warnings("the number of sheets in the input is different from PloGO2Results.xlsx")

	# merge plopathway results and input file
	list.mg = lapply(seq_len(length(sheets)), function(ii) 
		merge(readWorkbook(file.path(outFolder,"PloGO2Results.xlsx"), ii), data.list[[ii]], by.x="Protein", by.y=colName, 
		all=TRUE, sort=FALSE))

	wb <- createWorkbook()

	for(ii in seq_along(list.mg)) {
		addWorksheet(wb, sheets[ii])
		writeData(wb, sheets[ii], list.mg[[ii]])
		
	}

	addWorksheet(wb, sheets.plogo[length(sheets.plogo)])
	writeData(wb, sheets.plogo[length(sheets.plogo)],
		readWorkbook(file.path(outFolder,"PloGO2Results.xlsx"), length(sheets.plogo)) )

	saveWorkbook(wb, file.path(outFolder,"PloGO2Results.xlsx"), overwrite=TRUE)

	res
}
				