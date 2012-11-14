\name{DivideByPiR2}
\alias{DivideByPiR2}
\title{
	Divides a value by pi R^2
}
\description{
  This function is designed to simplify plotting K and similar functions dividing their values by  pi R^2
}
\usage{
DivideByPiR2(Value, r)
}
\arguments{
  \item{Value}{
  A vector of values.
 }
  \item{r}{
  The vector of distances
 }
}
\value{
  A vector of values.
}
\details{
  The classical normalization (Besag, 1977) of K(r) is L(r)=sqrt(K/pi). The normalization by pi R^2 is advocated by Marcon and Puech (2012) to allow the unification of distance-based measures of spatial structures.
}
\references{
 	Besag, J. E. (1977). Comments on Ripley's paper. Journal of the Royal Statistical Society B 39(2): 193-195.\cr
	Marcon, E. and F. Puech (2012). A typology of distance-based measures of spatial concentration. HAL SHS. 00679993.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
	\code{\link{K.r}}, \code{\link{L.r}}, \code{\link{Kmm.r}}, , \code{\link{D.r}}
}
\examples{
data(paracou16)
# Keep only 20% of points to run this example
X <- rthin(paracou16, 0.2)
plot(X)

# Calculate K
r <- 0:30
ActualValues.X <- K.r(X, r)

# Calculate confidence envelope (should be 1000 simulations, reduced to 4 to save time)
NumberOfSimulations <- 4
Alpha <- .10
LocalEnvelope.X <- KEnvelope(NumberOfSimulations, Alpha, X, r)

# Plot
PlotResults(r, DivideByPiR2(ActualValues.X, r), lapply(LocalEnvelope.X, DivideByPiR2, r), ylab="K / (pi R^2)", ReferenceValue=1)

}