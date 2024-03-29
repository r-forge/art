\name{powerTransform}
\alias{powerTransform}
\alias{powerTransform.default}
\alias{powerTransform.lm}
\alias{powerTransform.formula}

\title{Finding univariate or multivariate power transformations}
\description{
  \code{powerTransform} computes members of families of transformations indexed by one
parameter, the Box-Cox power family, or the Yeo and Johnson (2000) family, or the
basic power family, interpreting zero power as logarithmic. 
The family can be modified to have Jacobian one, or not, except for the basic
power family.
}
\usage{
powerTransform(object,...)

\S3method{powerTransform}{default}(object,...)

\S3method{powerTransform}{lm}(object, ...)

\S3method{powerTransform}{formula}(object, data, subset, weights, na.action,
  ...)
}

\arguments{
  \item{object}{For powerTransformation, this can either be an object of type
  \sQuote{lm}, a formula or a matrix or vector; see below.}
  \item{data}{A data frame or environment, as in \sQuote{lm}}
  \item{subset}{Case indices to be used, as in \sQuote{lm}.}
  \item{weights}{Case weights, is in \sQuote{lm}.}
  \item{na.action}{Missing value action, as in \sQuote{lm}.}
  \item{...}{Additional arguments that are passed to \code{\link{estimateTransform}}
  that does the actual computing, or the \code{\link{optim}} 
  function that does the maximization.  See the documentation for these functions
  for the arguments that are permitted, including \code{family} for
  setting the power transformation family.}
}

\details{The function powerTransform is used to estimate normalizing transformations
of a univariate or a multivariate random variable.  For a univariate transformation,
a formula like \code{z~x1+x2+x3} will find estimate a transformation for the response
\code{z} from the family of transformations indexed by the parameter \code{lambda}
that makes the residuals from the regression of the transformed \code{z} on the predictors
as closed to normally distributed as possible.  This generalizes the Box and
Cox (1964) transformations to normality only by allowing for families other than the
power transformations used in that paper.  

For a formula like \code{cbind(y1,y2,y3)~x1+x2+x3}, the three variables on
the left-side are all transformed, generally with different transformations
to make all the residuals as close to
normally distributed as possible.  \code{cbind(y1,y2,y3)~1} would specify transformations
to multivariate normality with no predictors.  This generalizes the multivariate
power transformations suggested by Velilla (1993) by allowing for different
families of transformations, and by allowing for predictors.  Cook and Weisberg (1999) 
and Weisberg (2005) suggest the usefulness of transforming
a set of predictors \code{z1, z2, z3} for multivariate normality and for transforming
for multivariate normality conditional on levels of a factor, which is equivalent
to setting the predictors to be indicator variables for that factor.  

Specifying the first argument as a vector, for example
\code{powerTransform(ais$LBM)}, is equivalent to
\code{powerTransform(LBM~1, ais)}.  Similarly,
\code{powerTransform(cbind(ais$LBM, ais$SSF))}, where the first argument is a matrix
rather than a formula is equivalent to
\code{powerTransform(cbind(LBM,SSF)~1, ais)}. 

Two families of power transformations are available.
The bcPower family of \emph{scaled power transformations},
\code{family="bctrans"}, 
equals \eqn{(U^{\lambda}-1)/\lambda}{(U^(lambda)-1)/lambda} 
for \eqn{\lambda}{lambda} \eqn{\neq}{not equal to} zero, and 
\eqn{\log(U)}{log(U)} if \eqn{\lambda =0}{lambda = 0}.
  
If \code{family="yjPower"} then the Yeo-Johnson transformations are used.  
This is is Box-Cox transformation of \eqn{U+1} for nonnegative values, 
and of \eqn{|U|+1} with parameter \eqn{2-\lambda}{2-lambda} for \eqn{U}
negative.
  
Other families can be added by writing a function whose first argument is a
matrix or vector to be transformed, whose second argument is the value of the
transformation parameter.  The function must return modified transformations
so that the Jacobian of the transformation is equal to one; see Cook and
Weisberg
(1982).

The function \code{powerTransform} is a front-end for
\code{\link{estimateTransform}}.  

The function \code{\link{testTransform}} is used to obtain likelihood ratio
tests for
any specified value for the transformation parameters.  It is used by the
summary method for powerTransform objects.
  
}
\value{
The result of \code{powerTransform} is an object of class \code{powerTransform}
that gives the estimates of the the transformation parameters and related 
statistics.  The print method for the object will display the estimates only;
the summary method provides both the estimates, standard errors, marginal Wald
confidence intervals and relevant likelihood ratio tests.  

Several helper functions are available.  The \code{coef} method returns 
the estimated transformation parameters, while \code{coef(object,round=TRUE)}
will
return the transformations rounded to nearby convenient values within 1.96
standard errors of the mle.
The \code{vcov} function returns the estimated covariance matrix of the
estimated
transformation parameters.  A \code{print} method is used to
print the objects and \code{summary} to provide more information.  By default
the summary method calls \code{testTransform} and provides likelihood ratio type
tests that all transformation parameters equal one and that all transformation
parameters equal 0, for log transformations, and for a convenient rounded value
not far from the mle.  The function can be called
directly to test any other value for lambda.
}
\references{Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations.
\emph{Journal
of the Royal Statisistical Society, Series B}. 26 211-46.

Cook, R. D. and Weisberg, S. (1999).  \emph{Applied Regression Including
Computing
and Graphics}.  Wiley.

Fox, J. (2002).  \emph{R and S-Plus Companion to Applied Regression}.  Sage. 

Velilla, S. (1993).  A note on the multivariate Box-Cox transformation to 
normality.  \emph{Statistics and Probability Letters}, 17, 259-263.

Weisberg, S. (2005) \emph{Applied Linear Regression}, third edition. Wiley.

Yeo, I. and Johnson, R. (2000).  A new family of
power transformations to improve normality or symmetry.  
\emph{Biometrika}, 87, 954-959.
}
\author{ Sanford Weisberg, sandy@stat.umn.edu }

\seealso{\code{\link{estimateTransform}},
\code{\link{optim}}, \code{\link{bcPower}}, \code{\link{transform}}. }
\examples{
data(ais,package="alr3")
# Box Cox Method, univariate
summary(p1 <- powerTransform(LBM~Ht+Wt+RCC,ais))
# fit linear model with transformed response:
coef(p1,round=TRUE)
summary(m1 <- lm(bcPower(LBM,p1$roundlam)~Ht+Wt+RCC,ais))
# Multivariate Box Cox
summary(powerTransform(cbind(LBM,Ht,Wt,RCC)~1,ais))
# Multivariate transformation to normality within levels of 'Sex'
summary(a3<-powerTransform(cbind(LBM,Ht,Wt,RCC)~Sex,ais))
# test lambda = (0 0 0 -1)
testTransform(a3,c(0,0,0,-1))
# save the rounded transformed values, plot them with a separate
# color for males and females
transformedY <- bcPower(with(ais,cbind(LBM,Ht,Wt,RCC)),coef(a3,round=TRUE))
\dontrun{pairs(transformedY,col=c("red","blue")[ais$Sex+1])}
# could also add the transformed variables to the data set using usual
# power transformations and the 'transform' function in the base package.
ais1 <- transform(ais,logLBM=log(LBM),logWt=log(Wt),logRCC=log(RCC))
}
\keyword{ regression}% at least one, from doc/KEYWORDS
