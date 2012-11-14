\name{SimulateM}
\alias{SimulateM}
\title{
	Simulations of point patterns to obtain values of M under the null hypothesis
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
	A string describing the null hypothesis to simulate. The null hypothesis may be "RandomLocation": points are redistributed on the actual locations; "RandomLabeling": randomizes point types, keeping locations and weights unchanged; "PopulationIndependence": keeps reference points unchanged, randomizes other point locations.
  }
  \item{CaseControl}{
	If TRUE, the case-control version of M is computed. ReferenceType points are cases, NeighborType points are controls.
  }
}

\value{
	A vector containing simulated M(r) values.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
	\code{\link{MEnvelope}}
}
\keyword{ internal }
