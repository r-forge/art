\name{bootCase}
\alias{bootCase}
\alias{bootCase.default}
\alias{bootCase.lm}
\alias{bootCase.glm}
\alias{bootCase.nls}
\alias{nextBoot}
\alias{nextBoot.lm}
\alias{nextBoot.nls}

\title{ Case bootstrap for regression models }
\description{
This routine does a case bootstrap resampling for regression models.  It
returns a matrix of the estimated coefficients from each of the bootstrap
samples.
}
\usage{
\S3method{bootCase}{lm}(object, f=coef, B=999)

\S3method{bootCase}{lm}(object, f=coef, B=999)

\S3method{bootCase}{nls}(object, f=coef, B=999)

nextBoot(object, sample)

}
\arguments{
  \item{object}{ A regression object of type \code{lm}, \code{glm} or class
\code{nls}.  May work with other regression objects that support the \code{update}
method and has a \code{subset} argument.  See details below.}
  \item{f}{A function that will be applied to the updated regression object to
compute the statistics of interest.  The default is \code{coef}, to return to
regression coefficient estimates.}
  \item{B}{Number of bootstrap samples.}
  \item{sample}{A sample of size n with replacement of the integers from 
1 to n that defines a bootstrap sample.}
}
\details{
This routine does the case-bootstrap described in the reference below.  Begin
with a regression object.  For each of B bootstrap samples, sample the rows of the data
matrix with replacement, and recompute and save estimates.  For nls objects there may
be convergence problems in the bootstrap.  The routine will continue until
convergence is attained B times, or until there are 25 consecutive failures to converge. 
\code{nextBoot} is an internal function that will update a model correctly, depending
on the class of the model object.

This simple routine should return a result with any S3 regression object that
can be updated using the \code{update} function and has a \code{subset}
argument.  In addition this method of resampling will not be
appropriate when the original data are grouped, so it is OK in general for
linear regression, logistic regression in which the response is either zero or
one.  With bionomial responses, one would generally want to resample one
observation, not all the observations in m trials, so this function will give
wrong answer.  The function can be used with Poisson regression with Poisson
sampling, but it is probably wrong for contingency tables with multinomial
sampling.  It is OK proportional odds models without Frequencies set, but
inappropriate with Frequencies.
}
\value{
A matrix with B rows and rank(object) columns giving the bootstrap estimates. 
These can be summarized as needed using standard R tools.
}
\references{
J. Fox and S. Weisberg (2011).  \emph{Companion to Applied Regression}, second edition.
Thousand Oaks:  Sage.

S. Weisberg (2005). \emph{Applied Linear Regression}, third edition.  New York: 
Wiley, Chapters 4 and 11.}

\author{ Sanford Weisberg, \email{sandy@stat.umn.edu}.  The error checking was
written by
Lexin Li.}

\seealso{See Also  \code{\link{update}}}


\examples{
m1 <- lm(Fertility ~ ., swiss)
betahat <- coef(m1)
betahat.boot <- bootCase(m1,B=99) # 99 bootstrap samples--too small to be useful
summary(betahat.boot)  # default summary
# bootstrap standard errors 
apply(betahat.boot,2,sd)
# bootstrap 95% confidence intervals
cl <- function(x) quantile(x,c(.025,.975))
apply(betahat.boot,2,cl)

}
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{internal}% __ONLY ONE__ keyword per line
