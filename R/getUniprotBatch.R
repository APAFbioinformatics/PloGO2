
# use Uniprot web site and RCurl
getUniprotBatch = function(filters="accession", values=IDList, attributes=c("id", "protein+names", "go-id")) {


query = paste("https://www.uniprot.org/uniprot/?query=", 
		paste(paste(filters, ":",values, sep=""), collapse="+OR+"), 
		"&columns=", paste(attributes, collapse=","), "&format=tab", sep="")
# res = try(getURI(query, ssl.verifypeer = FALSE))  # uses RCurl
res = try(GET(query))   # uses httr

if (inherits(res, "try-error")) stop("Error downloading records from Uniprot");
n = length(attributes)

# httr
res = content(res)
res.table = strsplit(res, "\n")[[1]]
# res.table = strsplit(res, "\n")[[1]]  # RCurl response is different

res.list = apply(data.frame(res.table), 1, FUN=function(v){
				w = strsplit(t(v), "\t")[[1]]
				if (length(w) < n) w = c(w, rep(" ", n-length(w))); w })  # pad with empty spaces

res.final = t(res.list)[-1,]
colnames(res.final) = res.list[,1]

res.final
}

