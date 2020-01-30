extractPairs <- function(res.list, ag.list, pair=names(ag.list)[1:2], targetFile) {

row.idx <- grep(targetFile, rownames(ag.list[[1]]))
if(length(row.idx) == 0) Error("Target category name not found in the abundance data");

counts <- data.frame(res.list[[pair[1]]]$counts[,1,drop=FALSE], 
	res.list[[pair[2]]]$counts[,1, drop=FALSE])


abundance <- data.frame(t(ag.list[[pair[1]]][row.idx,]), 
	t(ag.list[[pair[2]]][row.idx,]))

colnames(abundance) <- pair
rownames(counts) <- rownames(abundance)

data.frame(counts, abundance)


}


