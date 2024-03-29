\name{deltaMethod}
\alias{deltaMethod}
\alias{deltaMethod.default}
\alias{deltaMethod.lm} 
\alias{deltaMethod.nls}
\alias{deltaMethod.multinom} 
\alias{deltaMethod.polr} 

\title{Estimate and standard error of a nonlinear function of estimated regression 
coefficients}
\description{
\code{deltaMethod} is a generic function that uses the delta method to get a 
first-order approximate 
standard error for a nonlinear function of a vector of random variables
with known or estimated covariance matrix.
}
\usage{
deltaMethod(object, ...)

\method{deltaMethod}{default}(object, g ,var, func=g, ...)
\method{deltaMethod}{lm} (object, g, var=vcov,parameterPrefix="b",...)
\method{deltaMethod}{nls}(object, g, var=vcov,...)
\method{deltaMethod}{multinom} (object, g, var=vcov,parameterPrefix="b",...)
\method{deltaMethod}{polr} (object, g, var=vcov,...)

}
\arguments{
  \item{object}{For the default method, \code{object} is a named vector of \code{p} 
  elements.  This means that the call \code{names(object)} would return a list
  of \code{p} character strings that are the names of the elements of
  \code{object}.  For the other methods, \code{object} is a 
  regression object for which \code{coef(object)} returns a vector of parameter
  estimates.}
  \item{g}{A quoted string that is the function of the parameter estimates to be
  evaluated; see the details below. }
  \item{var}{The (estimated) covariance matrix of the coefficient
  estimates.  For the default method, this argument is required.  For all
  other methods, this argument must either provide the estimated covariance
  matrix or a function that when applied
  to \code{object} returns a covariance matrix.  The default is 
  to use the function \code{vcov}.}
  \item{func}{A quoted string used to annotate output.  The default of
  \code{func = g} is usually appropriate.}
  \item{parameterPrefix}{Typically a single letter with default
  \code{"b"} giving the prefix of the
  names of the parameter names used in the argument \code{g}; see details.}
  \item{...}{Additional arguments; not currently used.}
}

\details{
Suppose \eqn{x} is a random vector of length \eqn{p} that is at least approximately 
normally distributed with mean \eqn{\beta} and estimated covariance matrix 
\eqn{C}.   Then any function \eqn{g(\beta)} of \eqn{\beta}, is estimated by
\eqn{g(x)}, which is in large samples normally distributed with mean \eqn{g(\beta)}
and estimated variance \eqn{h'Ch}, where \eqn{h} is the first derivative of
\eqn{g(\beta)} with respect to \eqn{\beta} evaluated at \eqn{x}.  This function
returns both \eqn{g(x)} and its standard error, the square root of the estimated
variance.

The default method requires that you provide \eqn{x} in the argument
\code{object}, \eqn{C} in the argument \code{var}, and a text expression
in argument \code{g} that when evaluated gives the function \eqn{g}.  

Since
the delta method is often applied to functions of regression parameter 
estimates, the argument \code{object} may be the name of a regression
object from which the vector \eqn{x} will be taken from \code{coef(object)},
and \eqn{C} will be taken from \code{vcov(object)} unless you provide
some other estimate of variance, for example, using a sandwich estimator.  Methods
have been provided for many common regression models.   

For regression models for which methods are not provided, you must extract
the named vector of coefficient estimates and and estimate of its covariance
matrix and then apply the default deltaMethod function.

In the argument \code{g} you must provide a quoted character string
that gives the function of interest, for example \code{g="b1/b2"}, where
\code{b1} and \code{b2} are names of two of the coefficient estimates. 

For \code{lm}, \code{glm} and \code{polr} objects, the parameters are named for the
corresponding regressors.  For example, in in
\code{m2 <- lm(Y~X1+X2)}, \code{names(coef(m2))} returns
\code{c("(Intercept)","X1","X2")}.  If you wanted the estimate and standard 
error for the ratio of the coefficient estimates for \code{X1} and the
\code{(Intercept)}, you could either set \code{g = "X1/(Intercept)"} or
you could set \code{g = "b0/b1"}.  In general, you
can refer to the parameters using \dQuote{b0}, \dQuote{b1},...,\dQuote{bp}, 
if \code{parameterPrefix} is
left at its default value of \dQuote{b}.  In this case, \dQuote{b0} is the 
intercept (if the model has no intercept, then the numbering of the 
parameters starts with 1, not 0), \dQuote{b1} is the first estimated parameter
after the intercept, and so on.  You can also use the names of the corresponding
regressors in place of the parameter names.  for \code{multinom} objects you 
can only use the parameter names starting with the prefix, not the names of
the corresponding regressors. 

For nonlinear regression objects of type nls, the call \code{coef(object)} 
returns the estimated
coefficient vectors with names corresponding to parameter names.  
For example,           
\code{m2 <- nls(y~theta/(1+gamma*x), start = list(theta=2,gamma=3))} will
have parameters named  \code{c("theta", "gamma")}.   
In many other familiar regression methods, such as lm and glm, the names of
the coefficient estimates are the corresponding variable names, not
parameter names. 

Special characters should be avoided in the names of the elements of \eqn{x} 
as these can cause problems, and the \code{D} function used to compute derivatives may 
get confused.  However, embedded spaces or \dQuote{:} are permitted.   
}

\value{
A data.frame with two components
named \code{Estimate} for the estimate, \code{SE} for its standard error.
The value of \code{g} is given as a row label.  
}

\seealso{First derivatives of \code{g} are computed using symbolic differentiation
by the function \code{\link{D}}.}

\references{
J. Fox and S. Weisberg (2011). \emph{R Companion to Applied Regression},
second edition, Sage.

S. Weisberg (2005), \emph{Applied
Linear Regression}, third edition, Wiley, Section 6.1.2}
\author{Sanford Weisberg, \email{sandy@stat.umn.edu}
and John Fox,  \email{jfox@mcmaster.ca}}

\examples{
# cakes is a data frame with response Y, predictors X1 X2 
data(cakes,package="alr3") 
m1 <- lm(Y~ X2 + I(X2^2), data = cakes) # quadratic polynomial
deltaMethod(m1, "-b1/(2*b2)") # X2 that maximizes the quadratic
# second order polynomial in two predictors:
m2 <- lm(Y ~ X1 + X2 + I(X1^2) + I(X2^2) + X1:X2, data=cakes)
# Find X1 to maximize Y when X2=350:
deltaMethod(m2,"(b1+b5*350)/(-2*b3)")
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ regression }
