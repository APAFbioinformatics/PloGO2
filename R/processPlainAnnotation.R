processPlainAnnotation <- function (files, IDlist=NULL, datafile = NULL, printFiles = FALSE, 
    format = "compact", datafile.ignore.cols = 1, term.names=NULL) 
{

# if empty id list, generate it as the set of unique ID's from all the files

if (is.null(IDlist)) IDlist <- setdiff(unique(unlist(sapply(files, FUN=function(x){read.annot.file(x)[,2]}))), c(""," "));



    annotation.list <- lapply(files, FUN = function(fname) {
        res <- processAnnotFile(fname, IDlist, datafile, format = format, 
            datafile.ignore.cols = datafile.ignore.cols, term.names=term.names)
        fout <- basename(fname)
        if (printFiles) 
            writeGOannot(res, fname = paste("Annot", fout), datafile = datafile)
        res
    })

    names(annotation.list) <- basename(files)
    annotation.list
}


