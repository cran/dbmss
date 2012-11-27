\name{L.r}
\alias{L.r}
\title{
  Estimation of the \emph{L} function
}
\description{
  Estimates the \emph{L} function
}
\usage{
L.r(X, r, ReferenceType = "", NeighborType = "")
}
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
	\emph{L} is the normalized version of \emph{K}: \eqn{L(r)=\sqrt{\frac{K}{\pi}}-r}.
}
\value{
	A vector containing \emph{L(r)} values.
}
\references{
	Besag, J. E. (1977). Comments on Ripley's paper. \emph{Journal of the Royal Statistical Society B} 39(2): 193-195.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\note{
	\emph{L} was originally defined as \eqn{L(r)=\sqrt{\frac{K}{\pi}}}. It has been used as \eqn{L(r)=\sqrt{\frac{K}{\pi}}-r} in a part of the literature because this normalization is easier to plot.
}

\seealso{
	\code{\link{K.r}}, \code{\link{LEnvelope}}
}
\examples{

data(paracou16)
plot(paracou16)

# Calculate L
r <- 0:30
ActualValues.X <- L.r(paracou16, r)

# Plot
PlotResults(r, ActualValues.X, ylab="L", ReferenceValue=0)

}