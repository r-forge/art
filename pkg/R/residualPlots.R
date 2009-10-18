residualPlots <- function(m, ...){UseMethod("residualPlots")}
residualPlots.lm <- function(m,vars=~.,
     layout=NULL,ask,main="Residual Plots", fitted=TRUE,plot=TRUE,...){
  mf <- attr(model.frame(m), "terms")
  vform <- update(formula(m),vars)
  if(any(is.na(match(all.vars(vform), all.vars(formula(m))))))
     stop("Only predictors in the formula can be plotted.")
  terms <- attr(mf, "term.labels") # this is a list
  vterms <- attr(terms(vform), "term.labels")
# drop interactions (order > 1)
  vterms <- setdiff(vterms, terms[attr(mf, "order") > 1])
# keep only terms that are numeric or integer or factors or poly
  good <- NULL
  for (term in vterms) if(
      inherits(m$model[[term]], "numeric") |
      inherits(m$model[[term]], "integer") |
      inherits(m$model[[term]], "factor") | 
      inherits(m$model[[term]], "poly")) good <- c(good,term)
  nt <- length(good) + fitted
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
  ans <- NULL     
  if(!is.null(good)){
  for (term in good){
     nr <- nr+1
     ans <- rbind(ans,residualPlot(m,term,plot=plot,...))
     row.names(ans)[nr] <- term
    } }
  # Tukey's test
  if (fitted == TRUE){
   ans <- rbind(ans,residualPlot(m,"tukey",plot=plot,...))
   row.names(ans)[nr+1] <- "Tukey test"
   ans[nr+1,2] <- 2*pnorm(abs(ans[nr+1,1]),lower.tail=FALSE)}
  mtext(side=3,outer=TRUE,main, cex=1.4)
  dimnames(ans)[[2]] <- c("Test stat", "Pr(>|t|)")
  ans}
  
residualPlots.glm <- function(m, ...) {
 invisible(residualPlots.lm(m,...))
 }

residualPlot <- function(m,variable="tukey",type="pearson",
    plot=TRUE,add.quadratic=TRUE,
    identify.points = "xy",
    labels = names(residuals(m)[!is.na(residuals(m))]),
    id.n = 3, cex.identify=0.75, 
    col = palette()[2], col.lines = col[1],
    xlab, ylab, pch = 1, lwd = 2, ...){
 curvature <- class(m)[1] == "lm"
 string.capitalize <- function(string) {
     paste(toupper(substring(string,1,1)),substring(string,2),sep="")}
 ylab <- paste(string.capitalize(type),"Residuals")
 col <- match(variable,names(m$model))
 if(is.na(col) && variable != "tukey")
   stop(paste(variable,"is not a term in the mean function"))
 horiz <- if(variable == "tukey") predict(m) else m$model[[col]]
 lab <- if(variable == "tukey") {"Fitted values"} else variable
 ans <-
   if(inherits(horiz,"poly")) {
       horiz <- horiz[,1]
       lab <- paste("Linear part of",lab)
       c(NA,NA)}
   else if (class(horiz) == "factor") c(NA,NA)
   else if (curvature == TRUE) residCurvTest(m,variable)
   else  c(NA,NA)
# ans <- if (class(horiz) != "factor")  else c(NA,NA)
 if(plot==TRUE){
  plot(horiz,residuals(m,type=type),xlab=lab,ylab=ylab,...)
  abline(h=0,lty=2)
  if(class(horiz) != "factor") {
    if(add.quadratic==TRUE & curvature==TRUE){
        new <- seq(min(horiz),max(horiz),length=200)
        lm2 <- lm(residuals(m,type=type)~poly(horiz,2))
        lines(new,predict(lm2,list(horiz=new)),lty=3,lwd=2)
        }}}
  if (!is.logical(identify.points))
    showExtremes(horiz,residuals(m,type=type), labels=labels, ids=identify.points, 
      id.n=id.n, cex.id=cex.identify,res=residuals(m,type=type)) else 
      if (identify.points) identify(horiz,residuals(m,type=type), labels, cex=cex.identify)
  ans}
 
# September 24, 2009  Curvature testing is ONLY for lm's!
residCurvTest <- function(m,variable) {
 if(variable == "tukey") tukeyNonaddTest(m) else {
  if(is.na(match(variable, attr(m$terms,"term.labels"))))
     stop(paste(variable,"is not a term in the mean function")) else {
     xsqres <- qr.resid(m$qr,model.frame(m)[[variable]]^2)
     r <- residuals(m, type="pearson")
     m1 <- lm(r ~ xsqres, weights=weights(m))
     df.correction <- sqrt((df.residual(m)-1) / df.residual(m1))
     test <- summary(m1)$coef[2,3] * df.correction
     c(Test=test, Pvalue=2 * pt(-abs(test),df.residual(m)-1))
     }}}
residCurvTest.glm <- function(m) c(NA,NA)  # the above function might be OK...
     
tukeyNonaddTest <- function(m){
 qr <- m$qr
 fitsq <- qr.resid(qr,predict(m,type="response")^2)
 r <- residuals(m,type="pearson")
 m1 <- lm(r~fitsq,weights=weights(m))
 df.correction <- sqrt((df.residual(m)-1)/df.residual(m1))
 tukey <- summary(m1)$coef[2,3] * df.correction
 c(Test=tukey,Pvalue=2*pnorm(-abs(tukey)))
 }

