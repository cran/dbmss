\name{KinhomEnvelope}
\alias{KinhomEnvelope}
\title{
	Estimation of the confidence envelope of the \emph{Kinhom} function under its null hypothesis
}
\description{
	  Simulates point patterns according to the null hypothesis and returns the envelope of \emph{Kinhom} according to the confidence level.
}
\usage{
KinhomEnvelope(NumberOfSimulations, Alpha, X, r, ReferenceType="", 
    SimulationType="RandomPosition", lambda=NULL)
}
\arguments{
  \item{NumberOfSimulations}{
	The number of simulations to run.
  }
  \item{Alpha}{
	The risk level.
  }
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
  \item{SimulationType}{
	A string describing the null hypothesis to simulate. The null hypothesis, may be 
	"\emph{RandomPosition}": points are drawn in an inhomogenous Poisson process (intensity is either \code{lambda} or estimated from \code{X}); 
	"\emph{RandomLocation}": points are redistributed accross actual locations.
  }
  \item{lambda}{
	An estimation of the point pattern density, obtained by the \code{\link{density.ppp}} function.
  }
}
\details{
	This envelope is local, that is to say it is computed separately at each distance. See Loosmore and Ford (2006) for a discussion.
}
\value{
A list:
	\item{Simulations}{A matrix containing the simulated values (each line is a simulation, each column a value of \emph{Kinhom(R)}.}
	\item{Min}{A vector: the lower bound of the envelope.}
	\item{Max}{A vector: the upper bound of the envelope.}
}
\references{
	Kenkel, N. C. (1988). Pattern of Self-Thinning in Jack Pine: Testing the Random Mortality Hypothesis. \emph{Ecology} 69(4): 1017-1024.\cr
	Loosmore, N. B. and Ford, E. D. (2006). Statistical inference using the G or K point pattern spatial statistics. \emph{Ecology} 87(8): 1925-1931.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
	\code{\link{Kinhom.r}}
}

\examples{

data(paracou16)
# Keep only 20% of points to run this example
X <- rthin(paracou16, 0.2)
plot(X)

# Density of all trees
lambda <- density.ppp(X, bw.diggle(X))
plot(lambda)
V.americana <- X[X$marks$PointType=="V. Americana"]
plot(V.americana, add=TRUE)

# Calculate Kinhom according to the density of all trees
r <- 0:30
ActualValues.X <- Kinhom.r(X, r, "V. Americana", lambda)

# Calculate confidence envelope (should be 1000 simulations, reduced to 4 to save time)
NumberOfSimulations <- 4
Alpha <- .10
LocalEnvelope.X <- KinhomEnvelope(NumberOfSimulations, Alpha, X, r, , 
    SimulationType="RandomPosition", lambda=lambda)

# Plot
PlotResults(r, DivideByPiR2(ActualValues.X, r), lapply(LocalEnvelope.X, DivideByPiR2, r), 
    ylab="Kinhom / (pi R^2)", ReferenceValue=1)

}
