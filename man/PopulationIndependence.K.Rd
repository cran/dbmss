\name{PopulationIndependence.K}
\alias{PopulationIndependence.K}
\title{
	Simulations of a point pattern according to the null hypothesis of population independence defined for \emph{K}
}
\description{
  Simulates of a point pattern according to the null hypothesis of population independence defined for \emph{K}.
}
\usage{
PopulationIndependence.K(X, ReferenceType, NeighborType)
}
\arguments{
  \item{X}{
  A point pattern (\code{\link{ppp.object}}), marks must be a dataframe with two columns:\cr
    \code{PointType}: labels, as factors.\cr
    \code{PointWeight}: weights.
  }
  \item{ReferenceType}{
	One of the point types.
  }
  \item{NeighborType}{
	One of the point types.
  }
}
\details{
  Reference points are kept unchanged, other point positions are randomized by \code{\link{rshift}}.
}
\value{
  A new point pattern (an object of class \code{\link{ppp.object}}).
}
\references{
  Goreaud, F. et Pelissier, R. (2003). Avoiding misinterpretation of biotic interactions with the intertype K12 fonction: population independence vs random labelling hypotheses. \emph{Journal of Vegetation Science} 14(5): 681-692.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
  \code{\link{RandomPosition.K}}, \code{\link{PopulationIndependence.M}}, \code{\link{RandomLabeling.M}}
}
\examples{

# Simulate a point pattern with two types
X <- rpoispp(50) 
PointType   <- sample(c("A", "B"), X$n, replace=TRUE)
PointWeight <- runif(X$n, min=1, max=10)
X$marks <- data.frame(PointType, PointWeight)

par(mfrow=c(1,2))
plot(X, main="Original pattern, Point Type")

# Randomize it
Y <- PopulationIndependence.K(X, "A", "B")
# Invert the order of columns in mark to plot the point type, not the point weight
Y$marks <- data.frame(Y$marks$PointType, Y$marks$PointWeight)
# Points of type "A" are unchanged, points of type "B" have been moved altogether
plot(Y, main="Randomized pattern, Point Type")

}
