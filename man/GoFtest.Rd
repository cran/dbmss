\name{GoFtest}
\alias{GoFtest}
\title{
	Goodness of Fit test between a distance based measure of spatial structure and simulations of its null hypothesis
}
\description{
  Calculates the risk to reject the null hypothesis erroneously, based on the distribution of the simulations.
}
\usage{
GoFtest(ActualValues, SimulatedValues, r)
}
\arguments{
  \item{ActualValues}{
	A vector containing the values of a function (K, M,...) at several distances.
  }
  \item{SimulatedValues}{
  A matrix containing the simulated values (each line is a simulation, each column a value of a function(R)
  }
  \item{r}{
	A vector of distances.
  }
}
\details{
  This test was introduced by Diggle(1983) and extensively developped by Loosmore and Ford (2006) for K, and applied to M by Marcon et al. (2012).
}
\value{
  	A p-value.
}
\references{
  Diggle, P. J. (1983). Statistical analysis of spatial point patterns. Academic Press, London. 148 p.\cr
  Loosmore, N. B. and  Ford, E. D. (2006). Statistical inference using the G or K point pattern spatial statistics. Ecology 87(8): 1925-1931.\cr
  Marcon, E., F. Puech, et al. (2012). Characterizing the relative spatial structure of point patterns. International Journal of Ecology 2012(Article ID 619281): 11.  
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\note{
  No support exists in the litterature to apply the GoF test to non-cumulative functions (g, Kd...).\cr
  \code{\link{Ktest}} is a much better test (it does not rely on simulations) but it is limited to the K function against CSR in a rectangle window.
}

\seealso{
	\code{\link{Ktest}}
}
\examples{

# Simulate a Matern (Neyman Scott) point pattern
nclust <- function(x0, y0, radius, n) {
  return(runifdisc(n, radius, centre=c(x0, y0)))
}
X <- rNeymanScott(20, 0.2, nclust, radius=0.3, n=10)
plot(X)

# Calculate K
r <- seq(0, 0.3, 0.01)
ActualValues.X <- K.r(X, r)

# Calculate confidence envelope (should be 1000 simulations, reduced to 50 to save time)
NumberOfSimulations <- 50
Alpha <- .10
LocalEnvelope.X <- KEnvelope(NumberOfSimulations, Alpha, X, r)

# Plot
PlotResults(r, DivideByPiR2(ActualValues.X, r), lapply(LocalEnvelope.X, DivideByPiR2, r), ylab="K / (pi R^2)", ReferenceValue=1)

# GoF test. Power is correct if enough simulations are run (say >1000).
paste("p-value =", GoFtest(ActualValues.X, LocalEnvelope.X$Simulations, r))

}