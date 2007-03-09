deltaMethod <- function(object,g,parameter.prefix="b",var=vcov,...)
  {UseMethod("deltaMethod")}

### lm, glm, and others with unnamed parameters:
deltaMethod.default<-function(object,g,parameter.prefix="b",var=vcov,...)
{
   if(!is.character(g)) stop("The argument 'g' must be a character string")
   g<-parse(text=g)
   para <- coef(object)
   names(para) <- if ("(Intercept)" %in% names(para))
    paste(parameter.prefix,0:(length(para)-1),sep="") else
    paste(parameter.prefix,1:length(para),sep="")
   V <- var(object,...)             # vcov is the default
   computeDeltaMethod(V,g,para)}
   
# nls has named parameters so parameter.prefix is ignored
deltaMethod.nls<-function(object, g, parameter.prefix=NULL, var=vcov, ...)
{
   if(!is.character(g)) stop("The argument 'g' must be a character string")
   computeDeltaMethod(var(object,...),parse(text=g),coef(object))}
   
deltaMethod.drc<-function(object,g,parameter.prefix="b",var=vcov,...)
{
   if(!is.character(g)) stop("The argument 'g' must be a character string")
   para <- coef(object)
   names(para) <- paste(parameter.prefix,1:length(para),sep="")
   computeDeltaMethod(var(object,...),g=parse(text=g),para)}
   
# computes g evaluated at the data, and t(g')Vt(g'), the estimated
# standard error
computeDeltaMethod <- function(Var,g,values){
   q <- length(values)
   for(i in 1:q) {assign(names(values)[i], values[i])}
   est<-eval(g) 
   names(est)<-NULL
# derivative of function g of parameters
   gd<-NULL
   for(i in 1:q) {gd<-c(gd, eval(D(g, names(values)[i])))}
# compute se
   se.est<-as.vector(sqrt(t(gd) %*% Var %*% gd))
# output
   ans<-list(estimate=est, se=se.est, func=g)
   class(ans) <- c("deltaMethod","list")
   ans 
}

print.deltaMethod <- function(x,...){
     cat("Functions of parameters:  ")
     print.default(x$func)
     cat("Estimate =", x$est, "with se =", x$se, "\n")}
