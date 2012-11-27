\name{SimulateM}
\alias{SimulateM}
\title{
	Simulations of point patterns to obtain values of \emph{M} under the null hypothesis
}
\description{
	This function is used internally by \code{\link{MEnvelope}} to fill the simulation matrix.
}
\usage{
SimulateM(X, r, ReferenceType, NeighborType, SimulationType = "RandomLocation", CaseControl = FALSE)
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
  \item{SimulationType}{
	A string describing the null hypothesis to simulate. The null hypothesis may be "\emph{RandomLocation}": points are redistributed on the actual locations; "\emph{RandomLabeling}": randomizes point types, keeping locations and weights unchanged; "\emph{PopulationIndependence}": keeps reference points unchanged, randomizes other point locations.
  }
  \item{CaseControl}{
	Logical; if \code{TRUE}, the case-control version of \emph{M} is computed. \emph{ReferenceType} points are cases, \emph{NeighborType} points are controls.
  }
}

\value{
	A vector containing simulated \emph{M(r)} values.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
	\code{\link{MEnvelope}}
}
\keyword{ internal }
