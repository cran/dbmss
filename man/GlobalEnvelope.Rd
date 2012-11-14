\name{GlobalEnvelope}
\alias{GlobalEnvelope}
\title{
	Estimation of the global confidence interval of a matrix of simulations
}
\description{
  Calculates the global confidence interval envelope sensu Duaranton and Overman (2005) according to simulations of the null hypothesis of a function.
}
\usage{
GlobalEnvelope(Simulations, Alpha)
}
\arguments{
  \item{Simulations}{
  A matrix containing the simulated values (each line is a simulation, each column a value of a function(R)
  }
  \item{Alpha}{
	The risk level.
  }
}
\details{
  The global enveloppe is calculated by iteration: the simulations reaching one of the upper or lower values at any distance are eliminated at each step. The process is repeated until (Alpha / Number of simulations) simulations are dropped. The remaining upper and lower bounds at all distances constitute the global envelope. Interpolation is used if the exact ratio cannot be reached.\cr
  It is an alternative to the classical way to compute local confidence intervals (see \url{KEnvelope} for example.
}
\value{
A list:
	\item{Simulations}{A matrix containing the simulated values (each line is a simulation, each column a value of g(R)}
	\item{Min}{A vector: the lower bound of the envelope}
	\item{Max}{A vector: the upper bound of the envelope}
}
\references{
  Duranton, G. and Overman, H. G. (2005). Testing for Localisation Using Micro-Geographic Data. Review of Economic Studies 72(4): 1077-1106.  
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
  	\code{\link{KdEnvelope}}
}
\examples{

data(paracou16)
# Keep only 20% of points to run this example
X <- rthin(paracou16, 0.2)

# Calculate Kd
r <- 0:100
ActualValues.X <- Kd.r(X, r, "Q. Rosea", "Q. Rosea")

# Calculate confidence envelope (should be 1000 simulations, reduced to 20 to save time)
NumberOfSimulations <- 20
Alpha <- .20
LocalEnvelope.X <- KdEnvelope(NumberOfSimulations, Alpha, X, r, "Q. Rosea", "Q. Rosea")
GlobalEnvelope.X <- GlobalEnvelope(LocalEnvelope.X$Simulations, Alpha)

# Plot
PlotResults(r, ActualValues.X, LocalEnvelope.X, GlobalEnvelope.X, ylab="Kd", Legend=TRUE, LegendItems=c("Kd", "Local CI", "Global CI"), LegendPosition="bottomright")

}