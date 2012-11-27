\name{K.r}
\alias{K.r}
\title{
  Estimation of the \emph{K} function
}
\description{
	Estimates the \emph{K} function
}
\usage{
K.r(X, r, ReferenceType = "", NeighborType = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  A point pattern (\code{\link{ppp.object}}), marks must be a dataframe with two columns:\cr
    \code{PointType}: labels, as factors.\cr
    \code{PointWeight}: weights.
  }
  \item{r}{
	A vector of distances.
  }
  \item{ReferenceType}{
	One of the point types. Default is all point types.
  }
  \item{NeighborType}{
	One of the point types. Default is all point types.
  }
}
\details{
	\emph{K} is a cumulative, topographic measure of a point pattern structure.
}
\value{
	A vector containing \emph{K(r)} values.
}
\references{
	Ripley, B. D. (1976). The Foundations of Stochastic Geometry. \emph{Annals of Probability} 4(6): 995-998.\cr
	Ripley, B. D. (1977). Modelling Spatial Patterns. \emph{Journal of the Royal Statistical Society B} 39(2): 172-212.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\note{
	The computation of \code{K.r} relies on spatstat functions \code{\link{Kest}} and \code{\link{Kcross}}.
}
\seealso{
	\code{\link{L.r}}, \code{\link{Ktest}}, \code{\link{KEnvelope}}, \code{\link{Kest}}, \code{\link{Kcross}}
}
\examples{

data(paracou16)
plot(paracou16)

# Calculate K
r <- 0:30
ActualValues.X <- K.r(paracou16, r)

# Plot
PlotResults(r, DivideByPiR2(ActualValues.X, r), ylab="K / (pi R^2)", ReferenceValue=1)

}
