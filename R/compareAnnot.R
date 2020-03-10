compareAnnot <- function(res.list, referenceName, removeZeros = FALSE,
	correction = TRUE)
{
	file.list <- unlist(lapply(res.list, FUN = function(x) { x$fname }))
	reference <- grep(tolower(referenceName), tolower(file.list))

	if (length(reference) == 1) {
		Counts <- sapply(res.list, FUN = function(v) { v$counts[, 1] })
		Protein.counts <- sapply(res.list, FUN = function(v) { v$N })

		rest <- setdiff(1:ncol(Counts), reference)
		z.idx <- rep(FALSE, nrow(Counts))

		if (removeZeros) {
			z.idx <- (apply(Counts, 1, FUN = sum) == 0)
			Counts <- Counts[!z.idx, ]
		}

		# p.val <- matrix(NA, nrow(Counts), ncol(Counts)-1)
		p.val <- matrix(NA, nrow(Counts), ncol(Counts))

		# each category
		for (i in 1:nrow(Counts)) {
			for (j in rest) {
				mat <- matrix(c(Counts[i, j],
						Counts[i, reference],
						Protein.counts[j] - Counts[i, j],
						Protein.counts[reference] - Counts[i, reference]),
					byrow = TRUE, nrow = 2)

				#p.val[i] <- chisq.test(mat)$p.value
				if (min(mat) >= 3)
					p.val[i, j] <- fisher.test(mat)$p.value
			}
		}

		# p.val <- p.val[!z.idx,]
		GOIDlist <- rownames(res.list[[1]]$counts)
# DP 2012-05-24
# modify to allow for random lists, not just GO

		if (length(grep("^GO:", GOIDlist)) == length(GOIDlist)) {
			goTerms <- sapply(GOIDlist, FUN = function(v) { Term(GOTERM[[v]]) } )
		} else {
			goTerms <- GOIDlist
		}

		rownames(p.val) <- goTerms[!z.idx]
		# colnames(p.val) <- gsub(".*\\/","", file.names[rest])
#		colnames(p.val) <- gsub(".*\\/", "", file.list)
# TK 2011-10-25
		colnames(p.val) <- basename(file.list)

		if (correction) {
			for (i in 1:ncol(p.val))
				p.val[, i] <- p.adjust(p.val[, i], method = "fdr")
		}
	} else {
		# reference name not found or duplicate
		p.val <- NULL
	}

	p.val
}
