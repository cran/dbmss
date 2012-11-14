\name{KdEnvelope}
\alias{KdEnvelope}
\title{
	Estimation of the confidence envelope of the Kd function under its null hypothesis
}
\description{
	  Simulates point patterns according to the null hypothesis and returns the envelope of Kd according to the confidence level.
}
\usage{
KdEnvelope(NumberOfSimulations, Alpha, X, r, ReferenceType, NeighborType, Weighted = FALSE, SimulationType = "RandomLocation")
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
	One of the point types.
  }
  \item{NeighborType}{
	One of the point types.
  }
  \item{Weighted}{
	if TRUE, estimates the Kemp function.
  }
  \item{SimulationType}{
	A string describing the null hypothesis to simulate. The null hypothesis may be "RandomLocation": points are redistributed on the actual locations; "RandomLabeling": randomizes point types, keeping locations and weights unchanged; "PopulationIndependence": keeps reference points unchanged, randomizes other point locations.
  }
}
\details{
	This envelope is local, that is to say it is computed separately at each distance. See Loosmore and Ford (2006) for a discussion.
}
\value{
A list:
	\item{Simulations}{A matrix containing the simulated values (each line is a simulation, each column a value of Kd(R)}
	\item{Min}{A vector: the lower bound of the envelope}
	\item{Max}{A vector: the upper bound of the envelope}
}
\references{
	Kenkel, N. C. (1988). Pattern of Self-Thinning in Jack Pine: Testing the Random Mortality Hypothesis. Ecology 69(4): 1017-1024.\cr
	Loosmore, N. B. and Ford, E. D. (2006). Statistical inference using the G or K point pattern spatial statistics. Ecology 87(8): 1925-1931.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
	\code{\link{Kd.r}}
}

\examples{

data(paracou16)
# Keep only 50% of points to run this example
X <- rthin(paracou16, 0.5)
plot(X)

# Calculate Kd
r <- 0:100
ActualValues.X <- Kd.r(X, r, "Q. Rosea", "Q. Rosea")

# Calculate confidence envelope (should be 1000 simulations, reduced to 20 to save time)
NumberOfSimulations <- 20
Alpha <- .20
LocalEnvelope.X <- KdEnvelope(NumberOfSimulations, Alpha, X, r, "Q. Rosea", "Q. Rosea")
GlobalEnvelope.X <- GlobalEnvelope(LocalEnvelope.X$Simulations, Alpha)

# Plot
PlotResults(r, ActualValues.X, LocalEnvelope.X, GlobalEnvelope.X, ylab="Kd", Legend=TRUE, LegendItems=c("M", "Local CI", "Global CI"), LegendPosition="bottomright")

}