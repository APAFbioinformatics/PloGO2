
plotAbundanceBar = function(mat.abundance, mat.counts, min.count=5) {
  
mat.abundance[mat.abundance==0] = NA	

Max = apply(mat.counts, 1, FUN=max)

keep.idx = (Max > min.count) #( Max > 3 )
Objects = colnames(mat.abundance)

obj = Objects[keep.idx]


mat.abundance = mat.abundance[, (tolower(colnames(mat.abundance)) %in% tolower(obj)), drop=FALSE]
N = ncol(mat.abundance)

par(mar=c(1,1,1,1))
bar.abundance <- barplot(as.matrix(t(mat.abundance)), las=2, col=rainbow(N), horiz=TRUE, cex.names=.6)
legend("right", fill=rainbow(N), legend=colnames(mat.abundance), cex=0.6,xpd=TRUE, inset=-.1)

bar.abundance

}