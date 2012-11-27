\name{SimulateKinhom}
\alias{SimulateKinhom}
\title{
	Simulations of point patterns to obtain values of \emph{Kinhom} under the null hypothesis
}
\description{
	This function is used internally by \code{\link{KinhomEnvelope}} to fill the simulation matrix.
}
\usage{
SimulateKinhom(X, r, ReferenceType, SimulationType="RandomPosition", lambda)
}
\arguments{
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
  \item{SimulationType}{
	A string describing the null hypothesis to simulate. The null hypothesis, may be 
	"\emph{RandomPosition}": points are drawn in an inhomogenous Poisson process (intensity is either \code{lambda} or estimated from \code{X}); 
	"\emph{RandomLocation}": points are redistributed accross actual locations.
  }
  \item{lambda}{
	An estimation of the point pattern density, obtained by the \code{\link{density.ppp}} function.
  }
}

\value{
	A vector containing simulated \emph{Kinhom(r)} values.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
	\code{\link{KinhomEnvelope}}
}
\keyword{ internal }
