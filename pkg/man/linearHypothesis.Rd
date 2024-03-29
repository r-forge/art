%-------------------------------------------------------------------------------
% Revision history:
% checked in 2008-12-29 by J. Fox (corresponds to version 1.2-10 of car
%		with function renamed from linear.hypothesis)
%	2009-01-16 updated doc to correspond to changes in linearHypothesis. J. Fox
%-------------------------------------------------------------------------------

\name{linearHypothesis}
\alias{linearHypothesis}
\alias{lht}
\alias{linearHypothesis.lm}
\alias{linearHypothesis.glm}
\alias{linearHypothesis.mlm}
\alias{linearHypothesis.default}
\alias{print.linearHypothesis.mlm}

\title{Test Linear Hypothesis}

\description{
Generic function for testing a linear hypothesis, and methods
for linear models, generalized linear models, multivariate linear
models, and other models that have methods for \code{coef} and \code{vcov}. 
}

\usage{
linearHypothesis(model, ...)

lht(model, ...)

\method{linearHypothesis}{default}(model, hypothesis.matrix, rhs=NULL, 
    test=c("Chisq", "F"), vcov.=NULL, singular.ok=FALSE, verbose=FALSE, ...)  

\method{linearHypothesis}{lm}(model, hypothesis.matrix, rhs=NULL,
    test=c("F", "Chisq"), vcov.=NULL, 
	white.adjust=c(FALSE, TRUE, "hc3", "hc0", "hc1", "hc2", "hc4"), 
	singular.ok=FALSE, ...)

\method{linearHypothesis}{glm}(model,  ...)

\method{linearHypothesis}{mlm}(model, hypothesis.matrix, rhs=NULL, SSPE, V,
    test, idata, icontrasts=c("contr.sum", "contr.poly"), idesign, iterms, 
    P=NULL, title="", verbose=FALSE, ...)
       
\method{print}{linearHypothesis.mlm}(x, SSP=TRUE, SSPE=SSP, 
    digits=getOption("digits"), ...) 
}

\arguments{
  \item{model}{fitted model object. The default method works for models
   for which the estimated parameters can be retrieved by \code{coef} and
   the corresponding estimated covariance matrix by \code{vcov}. See the
   \emph{Details} for more information.}
  \item{hypothesis.matrix}{matrix (or vector) giving linear combinations
    of coefficients by rows, or a character vector giving the hypothesis
    in symbolic form (see \emph{Details}).}
  \item{rhs}{right-hand-side vector for hypothesis, with as many entries as
    rows in the hypothesis matrix; can be omitted, in which case it defaults
    to a vector of zeroes. For a multivariate linear model, \code{rhs} is a
    matrix, defaulting to 0.}
  \item{singular.ok}{if \code{FALSE} (the default), a model with aliased
    coefficients produces an error; if \code{TRUE}, the aliased coefficients
    are ignored, and the hypothesis matrix should not have columns for them.}
  \item{idata}{an optional data frame giving a factor or factors defining the
    intra-subject model for multivariate repeated-measures data. See 
    \emph{Details} for an explanation of the intra-subject design and for
    further explanation of the other arguments relating to intra-subject factors.}
  \item{icontrasts}{names of contrast-generating functions to be applied by default
    to factors and ordered factors, respectively, in the within-subject
    ``data''; the contrasts must produce an intra-subject model 
    matrix in which different terms are orthogonal.}
  \item{idesign}{a one-sided model formula using the ``data'' in \code{idata} and
    specifying the intra-subject design.}
  \item{iterms}{the quoted name of a term, or a vector of quoted names of terms,
   in the intra-subject design to be tested.}
  \item{P}{transformation matrix to be applied to the repeated measures in
    multivariate repeated-measures data; if \code{NULL} \emph{and} no 
    intra-subject model is specified, no response-transformation is applied; if
    an intra-subject model is specified via the \code{idata}, \code{idesign},
    and (optionally) \code{icontrasts} arguments, then \code{P} is generated 
    automatically from the \code{iterms} argument.}
  \item{SSPE}{in \code{linearHypothesis} method for \code{mlm} objects: 
    optional error sum-of-squares-and-products matrix; if missing,
    it is computed from the model. In \code{print} method for 
    \code{linearHypothesis.mlm} objects: if \code{TRUE}, 
    print the sum-of-squares and cross-products matrix for error.}
  \item{test}{character string, \code{"F"} or \code{"Chisq"},
    specifying whether to compute the finite-sample
    F statistic (with approximate F distribution) or the large-sample
    Chi-squared statistic (with asymptotic Chi-squared distribution). For a
    multivariate linear model, the multivariate test statistic to report --- one of
    \code{"Pillai"}, \code{"Wilks"}, \code{"Hotelling-Lawley"}, or \code{"Roy"}, 
    with \code{"Pillai"} as the default.}
  \item{title}{an optional character string to label the output.}
  \item{V}{inverse of sum of squares and products of the model matrix; if missing
    it is computed from the model.}
  \item{vcov.}{a function for estimating the covariance matrix of the regression
    coefficients, e.g., \code{\link{hccm}}, or an estimated covariance matrix
    for \code{model}. See also \code{white.adjust}.}        
  \item{white.adjust}{logical or character. Convenience interface to \code{hccm}
    (instead of using the argument \code{vcov.}). Can be set either to a character value
    specifying the \code{type} argument of \code{\link{hccm}} or \code{TRUE},
    in which case \code{"hc3"} is used implicitly. The default is \code{FALSE}.}
  \item{verbose}{If \code{TRUE}, the hypothesis matrix, right-hand-side
    vector (or matrix), and estimated value of the hypothesis 
    are printed to standard output; if \code{FALSE} (the default),
    the hypothesis is only printed in symbolic form and the value of the hypothesis
    is not printed.}
  \item{x}{an object produced by \code{linearHypothesis.mlm}.}
  \item{SSP}{if \code{TRUE} (the default), print the sum-of-squares and 
    cross-products matrix for the hypothesis and the response-transformation matrix.}
  \item{digits}{minimum number of signficiant digits to print.}
  \item{...}{aruments to pass down.}
}

\details{
  Computes either a finite sample F statistic or asymptotic Chi-squared
  statistic for carrying out a Wald-test-based comparison between a model
  and a linearly restricted model. The default method will work with any
  model object for which the coefficient vector can be retrieved by 
  \code{coef} and the coefficient-covariance matrix by \code{vcov} (otherwise
  the argument \code{vcov.} has to be set explicitely). For computing the
  F statistic (but not the Chi-squared statistic) a \code{df.residual} 
  method needs to be available. If a \code{formula} method exists, it is
  used for pretty printing.

  The method for \code{"lm"} objects calls the default method, but it
  changes the default test to \code{"F"}, supports the convenience argument
  \code{white.adjust} (for backwards compatibility), and enhances the output
  by residual sums of squares. For \code{"glm"} objects just the default
  method is called (bypassing the \code{"lm"} method).

  The function \code{lht} also dispatches to \code{linearHypothesis}.

  The hypothesis matrix can be supplied as a numeric matrix (or vector), 
  the rows of which specify linear combinations of the model coefficients, 
  which are tested equal to the corresponding entries in the righ-hand-side 
  vector, which defaults to a vector of zeroes. 

  Alternatively, the
  hypothesis can be specified symbolically as a character vector with one
  or more elements, each of which gives either a linear combination of
  coefficients, or a linear equation in the coefficients (i.e., with both
  a left and right side separated by an equals sign). Components of a
  linear expression or linear equation can consist of numeric constants, or
  numeric constants multiplying coefficient names (in which case the number
  precedes the coefficient, and may be separated from it by spaces or an
  asterisk); constants of 1 or -1 may be omitted.  Spaces are always optional.
  Components are separated by plus or minus signs. See the examples below.
  
  A linear hypothesis for a multivariate linear model (i.e., an object of
  class \code{"mlm"}) can optionally include an intra-subject transformation matrix
  for a repeated-measures design.
  If the intra-subject transformation is absent (the default), the multivariate 
  test concerns all of the corresponding coefficients for the response variables. 
  There are two ways to specify the transformation matrix for the
  repeated meaures: 
  \enumerate{
    \item The transformation matrix can be specified directly via the \code{P} 
    argument.
    
    \item A data frame can be provided defining the repeated-measures factor or
    factors 
    via \code{idata}, with default contrasts given by the \code{icontrasts}
    argument. An intra-subject model-matrix is generated from the one-sided formula 
    specified by the \code{idesign} argument; columns of the model matrix 
    corresponding to different terms in the intra-subject model must be orthogonal 
    (as is insured by the default contrasts). Note that the contrasts given in
    \code{icontrasts} can be overridden by assigning specific contrasts to the
    factors in \code{idata}. 
    The repeated-measures transformation matrix consists of the
    columns of the intra-subject model matrix corresponding to the term or terms
    in \code{iterms}. In most instances, this will be the simpler approach, and
    indeed, most tests of interests can be generated automatically via the
    \code{\link{Anova}} function.
  }
}

\value{
  For a univariate model, an object of class \code{"anova"} 
  which contains the residual degrees of freedom
  in the model, the difference in degrees of freedom, Wald statistic
  (either \code{"F"} or \code{"Chisq"}) and corresponding p value.
  
  For a multivariate linear model, an object of class 
  \code{"linearHypothesis.mlm"}, which contains sums-of-squares-and-product
  matrices for the hypothesis and for error, degrees of freedom for the 
  hypothesis and error, and some other information.
  
  The returned object normally would be printed.
}



\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
  
  Hand, D. J., and Taylor, C. C. (1987)
  \emph{Multivariate Analysis of Variance and Repeated Measures: A Practical
  Approach for Behavioural Scientists.} Chapman and Hall.
  
  O'Brien, R. G., and Kaiser, M. K. (1985)
  MANOVA method for analyzing repeated measures designs: An extensive primer.
  \emph{Psychological Bulletin} \bold{97}, 316--333.
}

\author{Achim Zeileis and John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[stats]{anova}}, \code{\link{Anova}}, \code{\link[lmtest]{waldtest}},
  \code{\link{hccm}}, \code{\link[sandwich]{vcovHC}}, \code{\link[sandwich]{vcovHAC}}, 
  \code{\link[stats]{coef}}, \code{\link[stats]{vcov}}

}

\examples{
data(Davis, package="car")
mod.davis <- lm(weight ~ repwt, data=Davis)

## the following are equivalent:
linearHypothesis(mod.davis, diag(2), c(0,1))
linearHypothesis(mod.davis, c("(Intercept) = 0", "repwt = 1"))
linearHypothesis(mod.davis, c("(Intercept)", "repwt"), c(0,1))
linearHypothesis(mod.davis, c("(Intercept)", "repwt = 1"))

## use asymptotic Chi-squared statistic
linearHypothesis(mod.davis, c("(Intercept) = 0", "repwt = 1"), test = "Chisq")


## the following are equivalent:
  ## use HC3 standard errors via white.adjust option
linearHypothesis(mod.davis, c("(Intercept) = 0", "repwt = 1"), 
    white.adjust = TRUE)
  ## covariance matrix *function*
linearHypothesis(mod.davis, c("(Intercept) = 0", "repwt = 1"), vcov = hccm)
  ## covariance matrix *estimate*
linearHypothesis(mod.davis, c("(Intercept) = 0", "repwt = 1"), 
    vcov = hccm(mod.davis, type = "hc3"))

data(Duncan, package="car")
mod.duncan <- lm(prestige ~ income + education, data=Duncan)

## the following are all equivalent:
linearHypothesis(mod.duncan, "1*income - 1*education = 0")
linearHypothesis(mod.duncan, "income = education")
linearHypothesis(mod.duncan, "income - education")
linearHypothesis(mod.duncan, "1income - 1education = 0")
linearHypothesis(mod.duncan, "0 = 1*income - 1*education")
linearHypothesis(mod.duncan, "income-education=0")
linearHypothesis(mod.duncan, "1*income - 1*education + 1 = 1")
linearHypothesis(mod.duncan, "2income = 2*education")

mod.duncan.2 <- lm(prestige ~ type*(income + education), data=Duncan)
coefs <- names(coef(mod.duncan.2))

## test against the null model (i.e., only the intercept is not set to 0)
linearHypothesis(mod.duncan.2, coefs[-1]) 

## test all interaction coefficients equal to 0
linearHypothesis(mod.duncan.2, coefs[grep(":", coefs)], verbose=TRUE)

## a multivariate linear model for repeated-measures data
## see ?OBrienKaiser for a description of the data set used in this example.

data(OBrienKaiser, package="car")

mod.ok <- lm(cbind(pre.1, pre.2, pre.3, pre.4, pre.5, 
                     post.1, post.2, post.3, post.4, post.5, 
                     fup.1, fup.2, fup.3, fup.4, fup.5) ~  treatment*gender, 
                data=OBrienKaiser)
coef(mod.ok)

## specify the model for the repeated measures:
phase <- factor(rep(c("pretest", "posttest", "followup"), c(5, 5, 5)),
    levels=c("pretest", "posttest", "followup"))
hour <- ordered(rep(1:5, 3))
idata <- data.frame(phase, hour)
idata
 
## test the four-way interaction among the between-subject factors 
## treatment and gender, and the intra-subject factors 
## phase and hour              
    
linearHypothesis(mod.ok, c("treatment1:gender1", "treatment2:gender1"),
    title="treatment:gender:phase:hour", idata=idata, idesign=~phase*hour, 
    iterms="phase:hour")
}

\keyword{htest}
\keyword{models}
\keyword{regression}

