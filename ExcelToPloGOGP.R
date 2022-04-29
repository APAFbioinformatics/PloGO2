############################################################
# parameters needed
# 	metafile name, by default "meta.csv" in local folder
#	metafile must contains files and groups
# defaults if any
#  	in and out folder by default local folder
#############################################################



ExcelToPloGOGP <- function(...) 
{

library(XLConnect)
library(openxlsx)
library(lattice)
library(PloGO)
# library(iTRAQRpac)

# initialize to defaults
termFile = NA
compareWithReference = "none"
data.file.name = "none"
colName = NA

args <- list(...)
for(i in 1:length(args)) {

flag <- substring(args[[i]], 0, 2)
value <- substring(args[[i]], 3, nchar(args[[i]]))

if(flag=='-f') fname <- value;
if(flag=='-c') colName <- value;
if(flag=='-t') termFile <- value;
if(flag=='-r') compareWithReference <- value;
if(flag=='-d') data.file.name <- value;

} 

res <- ExcelToPloGO( fname= fname, colName = colName,
	termFile = termFile, compareWithReference = compareWithReference, 
	data.file.name = data.file.name )
	
}

