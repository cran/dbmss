\name{RandomLabeling.M}
\alias{RandomLabeling.M}
\title{
	Simulations of a point pattern according to the null hypothesis of random labelling defined for \emph{M}
}
\description{
	Simulates of a point pattern according to the null hypothesis of random labelling defined for \emph{M}
}
\usage{
RandomLabeling.M(X)
}
\arguments{
  \item{X}{
  A point pattern (\code{\link{ppp.object}}), marks must be a dataframe with two columns:\cr
    \code{PointType}: labels, as factors.\cr
    \code{PointWeight}: weights.
  }
}
\details{
  Point types are randomized, keeping locations and weights unchanged.
}
\value{
  A new point pattern (an object of class \code{\link{ppp.object}}).
}
\references{
  Marcon, E. and Puech, F. (2010). Measures of the Geographic Concentration of Industries: Improving Distance-Based Methods. \emph{Journal of Economic Geography} 10(5): 745-762.\cr
  Marcon, E., F. Puech, et al. (2012). Characterizing the relative spatial structure of point patterns. \emph{International Journal of Ecology} 2012(Article ID 619281): 11.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
  \code{\link{PopulationIndependence.M}}, \code{\link{PopulationIndependence.K}}, \code{\link{RandomPosition.K}}
}
\examples{

# Simulate a point pattern with five types
X <- rpoispp(50) 
PointType   <- sample(c("A", "B", "C", "D", "E"), X$n, replace=TRUE)
PointWeight <- runif(X$n, min=1, max=10)
X$marks <- data.frame(PointType, PointWeight)

par(mfrow=c(2,2))
plot(X, main="Original pattern, Point Type")
X$marks <- data.frame(PointWeight=X$marks$PointWeight, PointType=X$marks$PointType)
plot(X, main="Original pattern, Point Weight")

# Randomize it
Y <- RandomLabeling.M(X)
Z <- Y
# Invert the order of columns in mark to plot the point type, not the point weight
Y$marks <- data.frame(PointType=Y$marks$PointType, PointWeight=Y$marks$PointWeight)
# Labels have been redistributed randomly accross locations
plot(Y, main="Randomized pattern, Point Type")
# But weights are unchanged
Y <- Z
plot(Y, main="Randomized pattern, Point Weight")

}