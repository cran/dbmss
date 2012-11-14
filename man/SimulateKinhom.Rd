\name{SimulateKinhom}
\alias{SimulateKinhom}
\title{
	Simulations of point patterns to obtain values of Kinhom under the null hypothesis
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
	"RandomPosition": points are drawn in an inhomogenous Poisson process (intensity is either lambda or estimated from X); 
	"RandomLocation": points are redistributed accross actual locations.
  }
  \item{lambda}{
	An estimation of the point pattern density, obtained by the \code{\link{density.ppp}} function.
  }
}

\value{
	A vector containing simulated Kinhom(r) values.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
	\code{\link{KinhomEnvelope}}
}
\keyword{ internal }
