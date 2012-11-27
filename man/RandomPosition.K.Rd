\name{RandomPosition.K}
\alias{RandomPosition.K}
\title{
	Simulations of a point pattern according to the null hypothesis of random position defined for \emph{K}
}
\description{
	Simulations of a point pattern according to the null hypothesis of random position defined for \emph{K}.
}
\usage{
RandomPosition.K(X)
}
\arguments{
  \item{X}{
  A point pattern (\code{\link{ppp.object}}), marks must be a dataframe with two columns:\cr
    \code{PointType}: labels, as factors.\cr
    \code{PointWeight}: weights.
  }
}
\details{
  Points marks are kept unchanged and their position is drawn in a binomial process by \code{\link{runifpoint}}.
}
\value{
  A new point pattern (an object of class \code{\link{ppp.object}}).
}
\note{
  Simulations in a binomial process keeps the same number of points, so that marks can be redistributed. If a real CSR simulation is needed and marks are useless, use \code{\link{rpoispp}}
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
  \code{\link{PopulationIndependence.K}}, \code{\link{PopulationIndependence.M}}, \code{\link{RandomLabeling.M}}
}
\examples{

# Simulate a point pattern with two types
X <- rpoispp(5) 
PointType   <- sample(c("A", "B"), X$n, replace=TRUE)
PointWeight <- runif(X$n, min=1, max=10)
X$marks <- data.frame(PointType, PointWeight)

par(mfrow=c(1,2))
plot(X, main="Original pattern, Point Type")

# Randomize it
Y <- RandomPosition.K(X)
# Invert the order of columns in mark to plot the point type, not the point weight
Y$marks <- data.frame(Y$marks$PointType, Y$marks$PointWeight)
# Points are randomly distributed
plot(Y, main="Randomized pattern, Point Type")

}