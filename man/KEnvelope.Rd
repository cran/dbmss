\name{KEnvelope}
\alias{KEnvelope}
\title{
	Estimation of the confidence envelope of the \emph{K} function under its null hypothesis
}
\description{
	  Simulates point patterns according to the null hypothesis and returns the envelope of \emph{K} according to the confidence level.
}
\usage{
KEnvelope(NumberOfSimulations, Alpha, X, r, ReferenceType = "", NeighborType = "", 
    SimulationType = "RandomPosition")
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
  \item{NeighborType}{
	One of the point types. Default is all point types.
  }
  \item{SimulationType}{
	A string describing the null hypothesis to simulate. The null hypothesis may be "\emph{RandomPosition}": points are drawn in a Poisson process; "\emph{RandomLabeling}": randomizes point types, keeping locations unchanged; "\emph{PopulationIndependence}": keeps reference points unchanged, randomizes other point locations.
  }
}
\details{
	This envelope is local, that is to say it is computed separately at each distance. See Loosmore and Ford (2006) for a discussion.
}
\value{
A list:
	\item{Simulations}{A matrix containing the simulated values (each line is a simulation, each column a value of \emph{K(R)}}
	\item{Min}{A vector: the lower bound of the envelope}
	\item{Max}{A vector: the upper bound of the envelope}
}
\references{
	Kenkel, N. C. (1988). Pattern of Self-Thinning in Jack Pine: Testing the Random Mortality Hypothesis. \emph{Ecology} 69(4): 1017-1024.\cr
	Loosmore, N. B. and Ford, E. D. (2006). Statistical inference using the G or K point pattern spatial statistics. \emph{Ecology} 87(8): 1925-1931.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
	\code{\link{K.r}}
}

\examples{

data(paracou16)
# Keep only 20% of points to run this example
X <- rthin(paracou16, 0.2)
plot(X)

# Calculate K
r <- 0:30
ActualValues.X <- K.r(X, r)

# Calculate confidence envelope (should be 1000 simulations, reduced to 20 to save time)
NumberOfSimulations <- 20
Alpha <- .10
LocalEnvelope.X <- KEnvelope(NumberOfSimulations, Alpha, X, r)

# Plot
PlotResults(r, DivideByPiR2(ActualValues.X, r), lapply(LocalEnvelope.X, DivideByPiR2, r), 
    ylab="K / (pi R^2)", ReferenceValue=1)

}
