\name{PopulationIndependence.M}
\alias{PopulationIndependence.M}
\title{
	Simulations of a point pattern according to the null hypothesis of population independence defined for M
}
\description{
  	Simulates of a point pattern according to the null hypothesis of population independence defined for M
}
\usage{
PopulationIndependence.M(X, ReferenceType)
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
}
\details{
  Reference points are kept unchanged, other point are redistributed randomly accross locations.
}
\value{
  A new point pattern (an object of class \code{\link{ppp.object}}).
}
\references{
  Marcon, E. and Puech, F. (2010). Measures of the Geographic Concentration of Industries: Improving Distance-Based Methods. Journal of Economic Geography 10(5): 745-762.\cr
  Marcon, E., F. Puech, et al. (2012). Characterizing the relative spatial structure of point patterns. International Journal of Ecology 2012(Article ID 619281): 11.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
  \code{\link{RandomLabeling.M}}, \code{\link{PopulationIndependence.K}}, \code{\link{RandomPosition.K}}
}
\examples{

# Simulate a point pattern with five types
X <- rpoispp(50) 
PointType   <- sample(c("A", "B", "C", "D", "E"), X$n, replace=TRUE)
PointWeight <- runif(X$n, min=1, max=10)
X$marks <- data.frame(PointType, PointWeight)

par(mfrow=c(1,2))
plot(X, main="Original pattern, Point Type")

# Randomize it
Y <- PopulationIndependence.M(X, "A")
# Invert the order of columns in mark to plot the point type, not the point weight
Y$marks <- data.frame(Y$marks$PointType, Y$marks$PointWeight)
# Points of type "A" (circles) are unchanged, all other points have been redistributed randomly accross locations
plot(Y, main="Randomized pattern, Point Type")

}