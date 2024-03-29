\name{leveneTest}
\alias{leveneTest}
\alias{leveneTest.formula}
\alias{leveneTest.lm}
\alias{leveneTest.default}

\title{Levene's Test}
\description{
 Computes Levene's test for homogeneity of variance across groups.
}
\usage{
leveneTest(y, ...)
\method{leveneTest}{formula}(y, data, ...)
\method{leveneTest}{lm}(y, ...)
\method{leveneTest}{default}(y, group, center=median, ...)
}

\arguments{
  \item{y}{response variable for the default method, or \code{lm} or 
  \code{formula} object. If \code{y} is a linear-model object or a formula,
  the variables on the right-hand-side of the model must all be factors and
  must be completely crossed.}
  \item{group}{factor defining groups.}
  \item{center}{The name of a function to compute the center of each group;
  	\code{mean} gives the original Levene's test; the default, \code{median}
  	provides a more robust test.}
  \item{data}{a data frame for evaluating the \code{formula}.}
  \item{...}{arguments to be passed down, e.g., \code{data} for the
  	\code{formula} and \code{lm} methods; can also be used to pass arguments to
  	the function given by \code{center} (e.g., \code{center=mean} and
  	\code{trim=0.1} specify the 10\% trimmed mean).}
}

\details{
  The version of Levene's test implemented here is based on absolute
  deviations from the group medians. Other, less robust, versions, not implemented,
  are based on absolute deviations from group means or trimmed means.
}

\value{
    returns an object meant to be printed showing the results of the test.
}

\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}
     
\author{John Fox \email{jfox@mcmaster.ca}; original generic version 
    contributed by Derek Ogle}
\note{adapted from a response posted by Brian Ripley to the r-help email list.}

\examples{
data(Moore, package="car")
with(Moore, leveneTest(conformity, fcategory))
with(Moore, leveneTest(conformity, interaction(fcategory, partner.status)))
leveneTest(conformity ~ fcategory*partner.status, data=Moore)
leveneTest(lm(conformity ~ fcategory*partner.status, data=Moore))
leveneTest(conformity ~ fcategory*partner.status, data=Moore, center=mean)
leveneTest(conformity ~ fcategory*partner.status, data=Moore, center=mean, trim=0.1)
}        

\keyword{htest}
