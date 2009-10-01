# October 1, 2009  avPlots by S. Weisberg.  avPlot by John Fox
avPlots <- function(model,vars=~., layout=NULL, ask, main="Added Variable Plot", ...){
  vars <- if(is.character(vars)) paste("~",vars) else vars
  mf <- attr(model.frame(model), "terms")
  vform <- update(formula(model),vars)
  if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
     stop("Only predictors in the formula can be plotted.")
  terms <- attr(mf, "term.labels") # this is a list
  vterms <- attr(terms(vform), "term.labels")
# drop interactions (order > 1)
  vterms <- setdiff(vterms, terms[attr(mf, "order") > 1])
# keep only terms that are numeric or integer or factors or poly
  good <- NULL
  for (term in vterms) if(
      inherits(model$model[[term]], "numeric") |
      inherits(model$model[[term]], "integer") ) good <- c(good,term)
  nt <- length(good)
  if (nt == 0) stop("No plots specified, or all predictors are factors")
  if(is.null(layout)){
   layout <- switch(min(nt,9),c(1,1),c(1,2),c(2,2),c(2,2),c(3,2),c(3,2),
                                c(3,3),c(3,3),c(3,3))}
  nr <- 0
  op<-par(no.readonly=TRUE, oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
  ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
  on.exit(par(op))
  if(prod(layout) > 1)
    par(mfrow=layout,ask=ask)
        #mar=c(3.5,4,.5,2)+.1,mgp=c(2,1,0),
        #cex.lab=1.0,cex=0.7)# mai=c(.6,.6,.1,.1)}
    else par(mfrow=layout,ask=ask)
  for (term in good) avPlot(model, term,...)
  mtext(side=3,outer=TRUE,main, cex=1.4)
  invisible(0)
 }
 
avp <- function(...) avPlots(...)
 
avPlot <-  function(model, ...) UseMethod("avPlot")

avPlot.lm <-
function (model, variable,
    identify.points = "xy",
    labels = names(residuals(model)[!is.na(residuals(model))]),
    id.n = 3, cex.identify=0.75, 
    col = palette()[2], col.lines = col[1],
    xlab, ylab, pch = 1, lwd = 2,  ...)
{
    variable <- if (is.character(variable) & 1 == length(variable))
        variable
    else deparse(substitute(variable))
    mod.mat <- model.matrix(model)
    var.names <- colnames(mod.mat)
    var <- which(variable == var.names)
    if (0 == length(var))
        stop(paste(variable, "is not a column of the model matrix."))
    response <- response(model)
    responseName <- responseName(model)
    if (is.null(weights(model)))
        wt <- rep(1, length(response))
    else wt <- weights(model)
    res <- lsfit(mod.mat[, -var], cbind(mod.mat[, var], response),
        wt = wt, intercept = FALSE)$residuals
    xlab <- if(missing(xlab)) paste(var.names[var], "| others") else xlab
    ylab <- if(missing(ylab)) paste(responseName, " | others")  else ylab
    plot(res[, 1], res[, 2], xlab = xlab,
        ylab = ylab,
        col = col, pch = pch, ...)
    abline(lsfit(res[, 1], res[, 2], wt = wt), col = col.lines, lwd = lwd)
    if (!is.logical(identify.points))
      showExtremes(res[, 1],res[, 2], labels=labels, ids=identify.points, 
          id.n=id.n, cex.id=cex.identify,res=residuals(model)) else 
      if (identify.points) identify(res[, 1], res[, 2], labels, cex=cex.identify)
}

avPlot.glm<-function(model, variable, 
    identify.points = "xy",
    labels = names(residuals(model)[!is.na(residuals(model))]),
    id.n = 3, cex.identify=0.75, 
    col = palette()[2], col.lines = col[1],
    xlab, ylab, pch = 1, lwd = 2,  type=c("Wang", "Weisberg"), ...){
    #last modified 20 Feb 2002 by J. Fox
    type<-match.arg(type)
    variable<-if (is.character(variable) & 1==length(variable)) variable
        else deparse(substitute(variable))
    mod.mat<-model.matrix(model)
    var.names<-colnames(mod.mat)
    var<-which(variable==var.names)
    if (0==length(var)) stop(paste(variable,"is not a column of the model matrix."))
    response<-response(model)
    responseName<-responseName(model)
    wt<-model$prior.weights
    mod<-glm(response~mod.mat[,-var]-1, weights=wt, family=family(model))
    res.y<-residuals(mod, type="pearson")
    wt<-if (type=="Wang") wt*model$weights else wt
    res.x<-lsfit(mod.mat[,-var], mod.mat[,var], wt=wt,    
        intercept=FALSE)$residuals
    xlab <- if(missing(xlab)) paste(var.names[var], "| others") else xlab
    ylab <- if(missing(ylab)) paste(responseName, " | others") else ylab
    plot(res.x, res.y, xlab=xlab, 
        ylab=ylab, col=col, pch=pch,...)
    abline(lsfit(res.x, res.y, wt=wt), col=col.lines, lwd=lwd)

    if (!is.logical(identify.points))
      showExtremes(res.x,res.y, labels=labels, ids=identify.points, 
          id.n=id.n, cex.id=cex.identify,res=residuals(model)) else
      if (identify.points) identify(res.x, res.y, labels, cex=cex.identify)
    }

  
  
