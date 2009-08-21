\name{invTranPlot}
\alias{invTranPlot}
\alias{invTranEstimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Choose a predictor transformation visually or numerically }
\description{
\code{invTranPlot}
draws a two-dimensional scatterplot of \eqn{Y}{Y} versus 
\eqn{X}{X}, along with the OLS
fit from the regression of \eqn{Y}{Y} on 
\eqn{(X^{\lambda}-1)/\lambda}{(Y^(lam)-1)/lam}.  \code{invTranEstimate}
find the nonlinear least squares estimate of \eqn{\lambda}{lambda} and its
standard error.
}
\usage{
invTranPlot(x,y,lambda=c(-1,0,1),lty=1:(1+length(lambda)),lwd=2,
        col.lines=rainbow(length(lambda)+1,start=0.7, end=0.1),xlab=deparse(substitute(x)),
        ylab=deparse(substitute(y)),family="bcpower",optimal=TRUE,
        key="topleft",...)

invTranEstimate(x,y,family="bcpower",confidence=0.95,...)
}

\arguments{
  \item{x}{The predictor variable }
  \item{y}{The response variable }
  \item{lambda}{The powers used in the plot.  The optimal power than minimizes
  the residual sum of squares is always added unless optimal is FALSE. }
  \item{family}{The transformation family to use, \code{"bcpower"},
  \code{"yjpower"}, or a user-defined family.}
  \item{confidence}{returns a profile likelihood confidence interval for the optimal 
transformation with this confidence level.  If FALSE, no interval is returned.}
  \item{optimal}{Include the optimal value of lambda?}
  \item{lty}{line types corresponding to the powers}
  \item{lwd}{the width of the plotted lines, defaults to 2 times the standard}
  \item{col.lines}{color of the fitted lines corresponding to the powers}
  \item{key}{The default is \code{"topleft"}, in which case a legend is added to
the top left corner of the plot; other choices include \code{"bottomright"}.
If \code{key} is a vector of two coordinates, the legend is drawn at the coordinates
specified. If \code{key} is FALSE, no key is drawn; if TRUE, you can place the
key interactively by clicking on the plot.}
  \item{xlab}{Label for the horizontal axis.}

  \item{ylab}{Label for the vertical axis.}

  \item{\dots}{ additional arguments passed to the plot method.}
}
\value{
\code{invTranPlot}
  returns a graph and a data.frame with \eqn{\lambda}{lam} in the 
  first column, and the residual sum of squares from the regression
  for that \eqn{\lambda}{lam} in the second column.

  \code{invTranEstimate} returns a list with elements \code{lambda} for the
  estimate, \code{se} for its standard error, and \code{RSS}, the minimum
  value of the residual sum of squares.  
}
\seealso{ \code{\link{inverseResponsePlot}},\code{\link{optimize}}}

\references{
Weisberg, S. (2005).  \emph{Applied Linear Regression}, third edition.
New York:  Wiley. }
\author{Sanford Weisberg, \email{sandy@stat.umn.edu} }

\examples{
data(baeskel, package = "alr3")
with(baeskel,invTranPlot(Sulfur,Tension))
with(baeskel,invTranEstimate(Sulfur,Tension))
}    
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{regression}