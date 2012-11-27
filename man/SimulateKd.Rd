\name{SimulateKd}
\alias{SimulateKd}
\title{
	Simulations of point patterns to obtain values of \emph{Kd} under the null hypothesis
}
\description{
	This function is used internally by \code{\link{KdEnvelope}} to fill the simulation matrix.
}
\usage{
SimulateKd(X, r, ReferenceType, NeighborType, Weighted = FALSE, SimulationType = "RandomLocation")
}
%- maybe also 'usage' for other objects documented here.
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
  \item{NeighborType}{
	One of the point types.
  }
  \item{Weighted}{
	Logical; if \code{TRUE}, estimates the \emph{Kemp} function.
  }
  \item{SimulationType}{
	A string describing the null hypothesis to simulate. The null hypothesis may be "\emph{RandomLocation}": points are redistributed on the actual locations; "\emph{RandomLabeling}": randomizes point types, keeping locations and weights unchanged; "\emph{PopulationIndependence}": keeps reference points unchanged, randomizes other point locations.
  }
}

\value{
	A vector containing simulated \emph{Kd(r)} values.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
	\code{\link{KdEnvelope}}
}
\keyword{internal }
