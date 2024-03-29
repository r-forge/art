#-------------------------------------------------------------------------------
# Revision history:
# checked in 2009-09-28 by J. Fox (moved/renamed from car)
#-------------------------------------------------------------------------------

# Bonferroni test for an outlier (J. Fox)

outlierTest <- function(model, ...){
	UseMethod("outlierTest")
}

outlierTest.lm <- function(model, cutoff=0.05, n.max=10, order=TRUE, labels=names(rstudent), ...){
	rstudent <- rstudent(model)
	labels <- if(is.null(labels)) seq(along=rstudent) else labels
	if (length(rstudent) != length(labels)) 
		stop("Number of labels does not correspond to number of residuals.")
	df <- df.residual(model) - 1
	rstudent <- rstudent[!is.na(rstudent)]
	n <- length(rstudent)
	p <- if (class(model)[1] == "glm")
			2*(pnorm(abs(rstudent), lower.tail=FALSE))
		else 2*(pt(abs(rstudent), df, lower.tail=FALSE))
	bp <- n*p 
	ord <- if (order) order(bp) else 1:n
	ord <- ord[bp[ord] <= cutoff]
	result <- if (length(ord) == 0){
			which <- which.max(rstudent)
			list(rstudent=rstudent[which], p=p[which], bonf.p=bp[which], signif=FALSE, cutoff=cutoff)
		}
		else {
			if (length(ord) > n.max) ord <- ord[1:n.max]
			result <- list(rstudent=rstudent[ord], p=p[ord], bonf.p=bp[ord], signif=TRUE, cutoff=cutoff)
		}
	class(result)<-"outlierTest"
	result
}

print.outlierTest<-function(x, digits=5, ...){
	if (!x$signif){
		cat("\nNo studentized residuals with Bonferonni p <", x$cutoff)
		cat("\nLargest |rstudent|:\n")
	}
	bp <- x$bonf
	bp[bp > 1] <- NA
	table <- data.frame(rstudent=x$rstudent, 
		"unadjusted p-value"=signif(x$p, digits), "Bonferonni p"=signif(bp, digits), 
		check.names=FALSE)
	rownames(table) <- names(x$rstudent)
	print(table)
	invisible(x)
}
