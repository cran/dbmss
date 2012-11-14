\name{Kinhom.r}
\alias{Kinhom.r}
\title{
  Estimation of the inhomogenous K function
}
\description{
	Estimates the Kinhom function
}
\usage{
Kinhom.r(X, r, ReferenceType = "", lambda=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  A point pattern (\code{\link{ppp.object}}), marks must be a dataframe with two columns:\cr
    \code{PointType}: labels, as factors.\cr
    \code{PointWeight}: weights (ignored in computation).
  }
  \item{r}{
	A vector of distances.
  }
  \item{ReferenceType}{
	One of the point types. Default is all point types.
  }
  \item{lambda}{
	An estimation of the point pattern density, obtained by the \code{\link{density.ppp}} function.
  }
}
\details{
	Kinhom is a cumulative, topographic measure of an inhomogenous point pattern structure.\cr
	By default, density estimation is performed at points by \code{\link{density.ppp}} using the optimal bandwith (\code{\link{bw.diggle}}). It can be calculated separately (see example), including at pixels if the point pattern is too large for the default estimation to succeed, and provided as the argument \code{lambda}:
	Arbia et al. (2012) for example use another point pattern as a reference to estimate density.\cr
	Bivariate Kinhom is not currently supported.
}
\value{
	A vector containing Kinhom(r) values.
}
\references{
	Baddeley, A. J., J. Moller, et al. (2000). Non- and semi-parametric estimation of interaction in inhomogeneous point patterns. Statistica Neerlandica 54(3): 329-350.\cr
	Arbia, G., G. Espa, et al. (2012). Clusters of firms in an inhomogeneous space: The high-tech industries in Milan. Economic Modelling 29(1): 3-11.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\note{
	The computation of \code{Kinhom.r()} relies on spatstat functions \code{\link{Kinhom}}, \code{\link{density.ppp}} and \code{\link{bw.diggle}}.
}
\seealso{
	\code{\link{KinhomEnvelope}}, \code{\link{Kinhom}}
}
\examples{

data(paracou16)

# Density of all trees
lambda <- density.ppp(paracou16, bw.diggle(paracou16))
plot(lambda)
V.americana <- paracou16[paracou16$marks$PointType=="V. Americana"]
plot(V.americana, add=TRUE)

# Calculate Kinhom according to the density of all trees
r <- 0:30
ActualValues.X <- Kinhom.r(paracou16, r, "V. Americana", lambda)

# Plot
PlotResults(r, DivideByPiR2(ActualValues.X, r), ylab="Kinhom / (pi R^2)", ReferenceValue=1)

}
