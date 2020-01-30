
countAndAbundance <- function(x, main="Annotation and abundance plot", legend.text=NULL,
	args.legend = list(x = "topright")) {

if (( ncol(x) %% 2)) Error("The annotation abundance matrix should have an even number of columns");

ord <- order(x[,1], decreasing = TRUE)
ord <- ord[c(rev(seq(from = 1, to = length(ord), by = 2)), seq(from = 2, to = length(ord), by = 2))]

x <- x[ord,]

annotation <- x[,1:(ncol(x)/2)]
abundance <- x[,(1+(ncol(x)/2)):ncol(x)]

if (is.null(legend.text)) { legend.text <- gsub("(.*)\\.(.*)$", "\\1", colnames(annotation)) }


par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
par(mar = c(5.1, 2.1, 0.1, 0.25))
mp <- barplot(-t(as.matrix(annotation)),
	space = c(0, 2), beside = TRUE, horiz = TRUE,
	xlab = "Count", axisnames = FALSE, axes = FALSE)
	at <- axTicks(1)
axis(1, at = at, labels = -at)

axis(4, at = colMeans(mp) + 2, labels = rownames(x), tick = FALSE, line = -0.5,
		hadj = 0.5, las = 2)
title(main = main, outer = TRUE)

par(mar = c(5.1, 0.25, 0.1, 2.1))
barplot(t(as.matrix(abundance)),
	space = c(0, 2), beside = TRUE, horiz = TRUE,
	xlab = "Abundance", axisnames = FALSE, axes = FALSE,
	legend.text = legend.text, args.legend = args.legend)

# Turns 0.00 into 0
	at <- axTicks(1)
	axis(1, at = at, labels = at)



}


