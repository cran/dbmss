\name{KmmEnvelope}
\alias{KmmEnvelope}
\title{
	Estimation of the confidence envelope of the \emph{Lmm} function under its null hypothesis
}
\description{
	  Simulates point patterns according to the null hypothesis and returns the envelope of \emph{Lmm} according to the confidence level.
}
\usage{
KmmEnvelope(NumberOfSimulations, Alpha, X, r, ReferenceType = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NumberOfSimulations}{
	The number of simulations to draw.
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
	One of the point types. Others are ignored. Default is all point types.
  }
}
\details{
	This envelope is local, that is to say it is computed separately at each distance. See Loosmore and Ford (2006) for a discussion.
}
\value{
A list:
	\item{Simulations}{A matrix containing the simulated values (each line is a simulation, each column a value of \emph{Kmm(R)}}
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
	\code{\link{Kmm.r}}
}
\examples{

data(paracou16)
# Keep only 20% of points to run this example
X <- rthin(paracou16, 0.2)
plot(X)

# Calculate Kmm
r <- seq(0, 30, 2)
ActualValues.X <- Kmm.r(X, r)

# Calculate confidence envelope (should be 1000 simulations, reduced to 10 to save time)
NumberOfSimulations <- 10
Alpha <- .10
LocalEnvelope.X <- KmmEnvelope(NumberOfSimulations, Alpha, X, r)

# Plot
PlotResults(r, DivideByPiR2(ActualValues.X, r), lapply(LocalEnvelope.X, DivideByPiR2, r), 
    ylab="Kmm / (pi R^2)", ReferenceValue=1)

}
