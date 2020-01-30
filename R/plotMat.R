
plotMat <- function(z, number, cor=FALSE, log=FALSE, main="Abundance Levelplot", ...) {


cn <- colorRampPalette(c("white", "blue"), space="rgb")

if (log) {
z.idx <- apply(log(z),2, FUN=sum) == 0
} else {
z.idx <- apply(z,2, FUN=sum) == 0
}

# remove those categories not present
z <- z[,!z.idx]
number <- number[!z.idx]

# remove those categories with one protein 
# n.idx <- number > 1
# z <- z[,n.idx]
# number <- number[n.idx] 

# sort in decreasing order of number

ord <- order(number)

z <- z[, ord]
number <- number[ord]



# DP 2012-05-24
# modify to allow for random lists, not just GO
if (length(grep("^GO:", colnames(z))) == length(colnames(z)) ) {
	nn <- sapply(colnames(z), FUN=function(x){Term(GOTERM[[x]])})
} else {
	nn <- colnames(z)
}
colnames(z) <- paste(nn, number)
if (log) z <- log(z)
if (cor) z <- cor(z);
plot(levelplot(z, col.regions=cn(100), scales=list(x=list(rot=45, cex=1.1), y=list(cex=1.2)), 
	xlab="", ylab="", main=main, ...))

TRUE

}

